package enzyme

import (
	"github.com/rmcl/restriction-enzymes/constants"
)

type RestrictionBatch struct {
	Enzymes []Enzyme
}

func NewRestrictionBatch(enzymes ...Enzyme) RestrictionBatch {
	return RestrictionBatch{
		Enzymes: enzymes,
	}
}

func (restrictionBatch *RestrictionBatch) Add(enzyme Enzyme) {
	restrictionBatch.Enzymes = append(restrictionBatch.Enzymes, enzyme)
}

/*
Find the next recognition site in a sequence after the provided offset.
that is cut by any of the enzymes in the batch.
*/
func (restrictionBatch *RestrictionBatch) GetNextRecognitionSite(
	sequence string,
	offset int,
	isCircular bool,
) (int, *Enzyme, constants.Strand) {

	nextSitePosition := -1
	nextSiteStrand := constants.Watson
	var nextEnzyme *Enzyme

	for _, enzyme := range restrictionBatch.Enzymes {
		enzymeSitePosition, _, strand := enzyme.GetNextRecognitionSite(sequence, offset, isCircular)
		if enzymeSitePosition < 0 {
			continue
		}

		if nextSitePosition < 0 || enzymeSitePosition < nextSitePosition {
			nextSitePosition = enzymeSitePosition
			nextSiteStrand = strand
			nextEnzyme = &enzyme
		}
	}

	return nextSitePosition, nextEnzyme, nextSiteStrand

}
