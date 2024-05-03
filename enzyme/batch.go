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

func NewRestrictionBatch(enzymes ...Enzyme) RestrictionBatch {
	return RestrictionBatch{
		Enzymes: enzymes,
	}
}

func (restrictionBatch *RestrictionBatch) Add(enzyme ...Enzyme) {
	restrictionBatch.Enzymes = append(restrictionBatch.Enzymes, enzyme...)
	restrictionBatch.combinedRegex = nil
}

func (restrictionBatch *RestrictionBatch) AddBatch(batch RestrictionBatch) {
	restrictionBatch.Enzymes = append(restrictionBatch.Enzymes, batch.Enzymes...)
	restrictionBatch.combinedRegex = nil
}

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

type EnzymeCutSiteMatch struct {
	Enzyme    *Enzyme
	Positions []int // Should we change this to cut site? Right now its the start of the recognition site
	//Strand    constants.Strand
}

func (enzymeBatch *RestrictionBatch) Search(
	sequence string,
	isCircular bool,
) (map[string]EnzymeCutSiteMatch, error) {
	enzymeCutSites := make(map[string]EnzymeCutSiteMatch)

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
			if matchRecord, ok := enzymeCutSites[result.Enzyme.Name]; ok {
				matchRecord.Positions = append(matchRecord.Positions, result.Position)
				enzymeCutSites[result.Enzyme.Name] = matchRecord
			} else {
				enzymeCutSites[result.Enzyme.Name] = EnzymeCutSiteMatch{
					Enzyme:    result.Enzyme,
					Positions: []int{result.Position},
				}
			}

			lastSitePosition = result.Position + 1
		}
	}

	return enzymeCutSites, nil
}

type RecognitionSiteResult struct {
	Position int
	Enzyme   *Enzyme
	Strand   constants.Strand
}

/*
Find the next recognition site in a sequence after the provided offset.
that is cut by any of the enzymes in the batch.
*/
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

	results := []RecognitionSiteResult{
		{
			Position: match[0] + offset,
			Enzyme:   enzyme,
			Strand:   strand,
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
			results = append(results, RecognitionSiteResult{
				Position: match[0] + offset,
				Enzyme:   additionalEnzyme,
				Strand:   strand,
			})
		}
	}

	return results
}
