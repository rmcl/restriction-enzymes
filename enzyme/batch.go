package enzyme

import (
	"regexp"
	"sort"
	"strings"

	"github.com/rmcl/restriction-enzymes/constants"
)

type RestrictionBatch struct {
	Enzymes         []Enzyme
	combinedRegex   *regexp.Regexp
	fwdSiteToEnzyme map[string][]*Enzyme
	revSiteToEnzyme map[string][]*Enzyme
	orderedEnzymes  []*Enzyme
}

// Create a new restriction batch with the given enzymes.
func NewRestrictionBatch(enzymes ...Enzyme) RestrictionBatch {
	return RestrictionBatch{
		Enzymes: enzymes,
	}
}

// Add one or more enzymes to the restriction batch.
func (restrictionBatch *RestrictionBatch) Add(enzyme ...Enzyme) {
	restrictionBatch.Enzymes = append(restrictionBatch.Enzymes, enzyme...)
	restrictionBatch.combinedRegex = nil
}

// Add all the enzymes from another restriction batch to this one.
func (restrictionBatch *RestrictionBatch) AddBatch(batch RestrictionBatch) {
	restrictionBatch.Enzymes = append(restrictionBatch.Enzymes, batch.Enzymes...)
	restrictionBatch.combinedRegex = nil
}

// Return a mapping of enzyme name to a list of sites in the sequence it cuts.
// This returns the position of the cut site on the watson strand to mirror
// Biopython's interface.
//
// If isCircular is true, the sequence is treated as circular and the
// function will return a recognition site, even if it spans the
// beginning and end of the sequence.
func (enzymeBatch *RestrictionBatch) Search(sequence string, isCircular bool) (map[string][]int, error) {
	watsonCutSitesByEnzyme := make(map[string][]int)

	lastSitePosition := 0
	for lastSitePosition >= 0 {

		results := enzymeBatch.GetNextRecognitionSite(
			sequence,
			lastSitePosition,
			isCircular,
		)

		if results == nil {
			break
		}

		for _, result := range results {
			enzyme := result.Enzyme
			if enzymeCutSites, ok := watsonCutSitesByEnzyme[result.Enzyme.Name]; ok {
				enzymeCutSites = append(enzymeCutSites, result.WatsonCutIndex)
				watsonCutSitesByEnzyme[result.Enzyme.Name] = enzymeCutSites
			} else {
				watsonCutSitesByEnzyme[enzyme.Name] = []int{result.WatsonCutIndex}
			}

			lastSitePosition = result.RecognitionSiteIndex + 1
		}
	}

	return watsonCutSitesByEnzyme, nil
}

// Get the next recognition site in the sequence after the offset
// for any enzyme in the batch.
//
// If isCircular is true, the sequence is treated as circular and the
// function will return the next recognition site, even if it spans the
// beginning and end of the sequence.
func (restrictionBatch *RestrictionBatch) GetNextRecognitionSite(
	sequence string,
	offset int,
	isCircular bool,
) []RecognitionSiteResult {

	if offset > len(sequence) {
		return nil
	}

	if restrictionBatch.Enzymes == nil || len(restrictionBatch.Enzymes) == 0 {
		return nil
	}

	if restrictionBatch.combinedRegex == nil {
		restrictionBatch.buildCombinedPattern()
	}

	sequence = strings.ToLower(sequence)
	searchSequence := sequence[offset:]
	originalSequenceLength := len(sequence)

	if isCircular {

		arbitrary_extension_length := 10

		if len(sequence) > arbitrary_extension_length {
			searchSequence += sequence[:arbitrary_extension_length]
		} else {
			// If the sequence is shorter then the recognition site, we
			// need to append the sequence to itself to avoid an index
			// out of range error.
			searchSequence += sequence
		}
	}

	// lower case the search sequence
	searchSequence = strings.ToLower(searchSequence)

	searchOffset := 0
	for {
		match := restrictionBatch.combinedRegex.FindStringIndex(searchSequence[searchOffset:])
		if match == nil {
			return nil
		}

		matchStart := searchOffset + match[0] + offset

		// If the match is beyond the original sequence length, we
		// should ignore it
		if matchStart >= originalSequenceLength {
			return nil
		}

		results := restrictionBatch.getResultsForMatchPosition(
			matchStart,
			sequence,
			isCircular,
		)
		if len(results) > 0 {
			return results
		}

		searchOffset += match[0] + 1
		if searchOffset >= len(searchSequence) {
			return nil
		}
	}
}

func (restrictionBatch *RestrictionBatch) getResultsForMatchPosition(
	matchStart int,
	sequence string,
	isCircular bool,
) []RecognitionSiteResult {
	results := []RecognitionSiteResult{}
	matchedEnzymeNames := map[string]bool{}
	for _, candidate := range restrictionBatch.orderedEnzymes {
		remainingSequence := sequence[matchStart:]
		if isCircular {
			if len(sequence) > len(candidate.Site) {
				remainingSequence += sequence[:len(candidate.Site)-1]
			} else {
				remainingSequence += sequence
			}
		}

		if forwardMatch := candidate.RegexpFor.FindStringIndex(remainingSequence); forwardMatch != nil && forwardMatch[0] == 0 {
			watsonCutIndex, crickCutIndex := candidate.GetCutSitePositions(matchStart, constants.Watson)
			results = append(results, RecognitionSiteResult{
				RecognitionSiteIndex: matchStart,
				Enzyme:               candidate,
				Strand:               constants.Watson,
				WatsonCutIndex:       watsonCutIndex,
				CrickCutIndex:        crickCutIndex,
			})
			matchedEnzymeNames[candidate.Name] = true
		}

		if reverseMatch := candidate.RegexpRev.FindStringIndex(remainingSequence); reverseMatch != nil && reverseMatch[0] == 0 {
			if matchedEnzymeNames[candidate.Name] {
				continue
			}

			watsonCutIndex, crickCutIndex := candidate.GetCutSitePositions(matchStart, constants.Crick)
			results = append(results, RecognitionSiteResult{
				RecognitionSiteIndex: matchStart,
				Enzyme:               candidate,
				Strand:               constants.Crick,
				WatsonCutIndex:       watsonCutIndex,
				CrickCutIndex:        crickCutIndex,
			})
		}
	}

	return results
}

// Build a combined regular expression pattern for all the enzymes in the batch.
// this private method is called when the combinedRegex is nil and is used to
// build the combined pattern for all the enzymes in the batch. This is much
// faster than searching for each enzyme individually.
func (restrictionBatch *RestrictionBatch) buildCombinedPattern() {
	restrictionBatch.fwdSiteToEnzyme = make(map[string][]*Enzyme)
	restrictionBatch.revSiteToEnzyme = make(map[string][]*Enzyme)
	restrictionBatch.orderedEnzymes = nil

	if len(restrictionBatch.Enzymes) == 0 {
		restrictionBatch.combinedRegex = nil
		return
	}

	enzymes := append([]Enzyme(nil), restrictionBatch.Enzymes...)
	sort.Slice(enzymes, func(i, j int) bool {
		if enzymes[i].Length != enzymes[j].Length {
			return enzymes[i].Length > enzymes[j].Length
		}

		if enzymes[i].Site != enzymes[j].Site {
			return enzymes[i].Site < enzymes[j].Site
		}

		return enzymes[i].Name < enzymes[j].Name
	})

	patternString := "("
	for i := range enzymes {
		enzyme := &enzymes[i]
		restrictionBatch.orderedEnzymes = append(restrictionBatch.orderedEnzymes, enzyme)

		fwdExpStr := strings.ToLower(enzyme.RegexpFor.String())
		revExpStr := strings.ToLower(enzyme.RegexpRev.String())

		patternString += "" + fwdExpStr + "|" + revExpStr + "|"

		// Unfortunately, we need to strip the first 4 characters from the
		// the regular expression pattern because it is "(?i)"
		fwdSite := strings.ToLower(enzyme.RegexpFor.String()[4:])
		revSite := strings.ToLower(enzyme.RegexpRev.String()[4:])

		if restrictionBatch.fwdSiteToEnzyme[fwdSite] == nil {
			restrictionBatch.fwdSiteToEnzyme[fwdSite] = []*Enzyme{}
			restrictionBatch.revSiteToEnzyme[revSite] = []*Enzyme{}
		}

		restrictionBatch.fwdSiteToEnzyme[fwdSite] = append(restrictionBatch.fwdSiteToEnzyme[fwdSite], enzyme)
		restrictionBatch.revSiteToEnzyme[revSite] = append(restrictionBatch.revSiteToEnzyme[revSite], enzyme)
	}
	patternString = patternString[:len(patternString)-1] + ")"
	restrictionBatch.combinedRegex = regexp.MustCompile(patternString)
}
