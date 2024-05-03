package enzyme

import (
	"regexp"
	"strings"

	"github.com/rmcl/restriction-enzymes/constants"
)

type RestrictionBatch struct {
	Enzymes         []Enzyme
	combinedRegex   *regexp.Regexp
	fwdSiteToEnzyme map[string]*Enzyme
	revSiteToEnzyme map[string]*Enzyme
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

	match := restrictionBatch.combinedRegex.FindStringIndex(searchSequence)
	if match == nil {
		return nil
	}

	// If the match is beyond the original sequence length, we
	// should ignore it
	if match[0] >= originalSequenceLength {
		return nil
	}

	matchedSite := searchSequence[match[0]:match[1]]

	strand := constants.Watson
	enzyme := restrictionBatch.fwdSiteToEnzyme[matchedSite]
	if enzyme == nil {
		strand = constants.Crick
		enzyme = restrictionBatch.revSiteToEnzyme[matchedSite]
	}

	if enzyme == nil {
		// This should never happen
		return nil
	}

	watsonCutIndex, crickCutIndex := enzyme.GetCutSitePositions(match[0]+offset, strand)
	results := []RecognitionSiteResult{
		{
			RecognitionSiteIndex: match[0] + offset,
			Enzyme:               enzyme,
			Strand:               strand,

			WatsonCutIndex: watsonCutIndex,
			CrickCutIndex:  crickCutIndex,
		},
	}

	// Confirm that there are not multiple enzymes that match the same site
	// This is possible if the recognition site is 1 bp shorter than the additional
	// enzymes recognition site
	for i := 1; i < 5; i++ {
		endIndex := match[1] + i
		if endIndex >= len(searchSequence) {
			break
		}

		additionalSite := searchSequence[match[0]:endIndex]
		additionalEnzyme := restrictionBatch.fwdSiteToEnzyme[additionalSite]
		if additionalEnzyme != nil {
			watsonCutIndex, crickCutIndex := enzyme.GetCutSitePositions(match[0]+offset, strand)
			results = append(results, RecognitionSiteResult{
				RecognitionSiteIndex: match[0] + offset,
				Enzyme:               additionalEnzyme,
				Strand:               strand,

				WatsonCutIndex: watsonCutIndex,
				CrickCutIndex:  crickCutIndex,
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
	restrictionBatch.fwdSiteToEnzyme = make(map[string]*Enzyme)
	restrictionBatch.revSiteToEnzyme = make(map[string]*Enzyme)

	patternString := "("
	for _, enzyme := range restrictionBatch.Enzymes {
		patternString += "" + enzyme.RegexpFor.String() + "|" + enzyme.RegexpRev.String() + "|"

		fwdSite := strings.ToLower(enzyme.RegexpFor.String()[4:])
		revSite := strings.ToLower(enzyme.RegexpRev.String()[4:])

		restrictionBatch.fwdSiteToEnzyme[fwdSite] = &enzyme
		restrictionBatch.revSiteToEnzyme[revSite] = &enzyme
	}
	patternString = patternString[:len(patternString)-1] + ")"
	restrictionBatch.combinedRegex = regexp.MustCompile(patternString)
}
