package enzyme

import (
	"fmt"
	"regexp"

	"github.com/rmcl/restriction-enzymes/constants"
)

type EnzymeNumberOfCuts int

const (
	OneCut      EnzymeNumberOfCuts = 1
	TwoCuts     EnzymeNumberOfCuts = 2
	UnknownCuts EnzymeNumberOfCuts = 0
)

type EnzymeCutType string

const (
	BluntEnd  EnzymeCutType = "blunt"
	StickyEnd EnzymeCutType = "sticky"
)

type Enzyme struct {
	Name string

	Site      string
	Length    int
	Substrate string

	RegexpFor *regexp.Regexp
	RegexpRev *regexp.Regexp

	CutType EnzymeCutType

	OverhangLength   int
	OverhangSequence string

	NumberOfCuts EnzymeNumberOfCuts

	FivePrimeCutSite  int
	ThreePrimeCutSite int

	FivePrimeCutSite2  int
	ThreePrimeCutSite2 int

	RebaseId                int
	InactivationTemperature int
	OptimalTemperature      int

	Uri string

	Suppliers []string

	References []string
}

/*
Find the next recognition site in a sequence after the provided offset.

If the enzyme is circular, the sequence is treated as circular and the
function will return the next recognition site, even if it spans the
beginning and end of the sequence.
*/
func (enzyme *Enzyme) GetNextRecognitionSite(
	sequence string,
	offset int,
	isCircular bool,
) []RecognitionSiteResult {

	remainingSequence := sequence[offset:]

	// If the enzyme is circular, we need to support the case where the
	// recognition site spans the beginning and end of the sequence.
	// To do this, we append the sequence to itself of length site - 1.
	if isCircular {
		if len(sequence) > len(enzyme.Site) {
			fmt.Println(enzyme)
			fmt.Println("sequence", sequence, len(sequence), len(enzyme.Site))
			remainingSequence += sequence[:len(enzyme.Site)-1]
		} else {
			// If the sequence is shorter then the recognition site, we
			// need to append the sequence to itself to avoid an index
			// out of range error.
			remainingSequence += sequence

		}
	}

	var watsonMatchIndex int
	var crickMatchIndex int

	watsonMatch := enzyme.RegexpFor.FindStringIndex(remainingSequence)
	if watsonMatch != nil {
		watsonMatchIndex = watsonMatch[0] + offset
	} else {
		watsonMatchIndex = -1
	}

	crickMatch := enzyme.RegexpRev.FindStringIndex(remainingSequence)
	if crickMatch != nil {
		crickMatchIndex = crickMatch[0] + offset
	} else {
		crickMatchIndex = -1
	}

	if watsonMatchIndex > -1 && crickMatchIndex > -1 {

		if watsonMatchIndex <= crickMatchIndex {
			return []RecognitionSiteResult{
				{
					Enzyme:   enzyme,
					Position: watsonMatchIndex,
					Strand:   constants.Watson,
				},
			}
		} else {
			return []RecognitionSiteResult{
				{
					Enzyme:   enzyme,
					Position: crickMatchIndex,
					Strand:   constants.Crick,
				},
			}
		}

	} else if watsonMatchIndex > -1 {
		return []RecognitionSiteResult{
			{
				Enzyme:   enzyme,
				Position: watsonMatchIndex,
				Strand:   constants.Watson,
			},
		}
	} else if crickMatchIndex > -1 {
		return []RecognitionSiteResult{
			{
				Enzyme:   enzyme,
				Position: crickMatchIndex,
				Strand:   constants.Crick,
			},
		}
	}

	// No match found
	return nil
}
