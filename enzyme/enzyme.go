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

	References []string
}

// A struct to hold a single recognition site result in a sequence.
type RecognitionSiteResult struct {
	Enzyme *Enzyme

	// The position of the start of the recognition site in the sequence.
	RecognitionSiteIndex int
	Strand               constants.Strand

	// The position of the cut site on the watson strand in the sequence.
	WatsonCutIndex int
	// The position of the cut site on the watson strand in the sequence.
	CrickCutIndex int
}

// Get the next recognition site in the sequence after the offset
//
// If isCircular is true, the sequence is treated as circular and the
// function will return the next recognition site, even if it spans the
// beginning and end of the sequence.
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

	var strand constants.Strand
	var position int

	if watsonMatchIndex > -1 && crickMatchIndex > -1 {
		if watsonMatchIndex <= crickMatchIndex {
			strand = constants.Watson
			position = watsonMatchIndex
		} else {
			strand = constants.Crick
			position = crickMatchIndex
		}
	} else if watsonMatchIndex > -1 {
		strand = constants.Watson
		position = watsonMatchIndex
	} else if crickMatchIndex > -1 {
		strand = constants.Crick
		position = crickMatchIndex
	} else {
		// No match found
		return nil
	}

	watsonIndex, crickIndex := enzyme.GetCutSitePositions(watsonMatchIndex, strand)
	return []RecognitionSiteResult{
		{
			Enzyme: enzyme,

			RecognitionSiteIndex: position,
			Strand:               strand,

			WatsonCutIndex: watsonIndex,
			CrickCutIndex:  crickIndex,
		},
	}
}

// Given a recognition site index and a strand, return the position of the cut sites
// on the watson and crick strands. The first int is the watson strand index and the second
// is the crick strand index.
func (enzyme *Enzyme) GetCutSitePositions(recognitionSiteIndex int, strand constants.Strand) (int, int) {
	var watsonCutIndex, crickCutIndex int

	if strand == constants.Watson {
		watsonCutIndex = recognitionSiteIndex + enzyme.FivePrimeCutSite
		crickCutIndex = recognitionSiteIndex + enzyme.ThreePrimeCutSite
	} else {
		watsonCutIndex = recognitionSiteIndex + enzyme.Length - enzyme.ThreePrimeCutSite
		crickCutIndex = recognitionSiteIndex + enzyme.Length - enzyme.FivePrimeCutSite
	}

	return watsonCutIndex, crickCutIndex
}

// Return a list of cut sites for the enzyme. This returns the position of the cut site
// on the watson strand to mirror Biopython
//
// If isCircular is true, the sequence is treated as circular and the function will
// return a recognition site even if it spans the beginning and end of the sequence.
func (enzyme *Enzyme) Search(
	sequence string,
	isCircular bool,
) ([]int, error) {

	watsonCutSites := make([]int, 0)

	lastSitePosition := 0
	for lastSitePosition >= 0 {

		results := enzyme.GetNextRecognitionSite(
			sequence,
			lastSitePosition,
			isCircular,
		)

		if results == nil {
			break
		}

		for _, result := range results {
			watsonCutSites = append(watsonCutSites, result.WatsonCutIndex)

			lastSitePosition = result.RecognitionSiteIndex + 1
		}
	}

	return watsonCutSites, nil
}
