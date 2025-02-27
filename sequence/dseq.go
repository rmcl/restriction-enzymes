/*
Sequence package provides a Dseq struct that represents a double stranded DNA sequence.
*/
package sequence

import (
	"fmt"

	"github.com/bebop/poly/transform"
	"github.com/rmcl/restriction-enzymes/constants"
	"github.com/rmcl/restriction-enzymes/enzyme"
)

type Dseq struct {
	// A string representing the watson (sense) DNA strand.
	Watson string

	// A string representing the crick (antisense) DNA strand.
	Crick string

	// A positive or negative number to describe the stagger between
	// start of the watson and crick strands.
	Overhang int

	// A string representing the geometry of the DNA sequence. Can be "linear" or "circular".
	Geometry constants.SequenceGeometry
}

func (dSeq *Dseq) Print() {
	if dSeq.Overhang > 0 {
		for i := 0; i < dSeq.Overhang; i++ {
			fmt.Print(" ")
		}
	}
	fmt.Println(dSeq.Watson)
	if dSeq.Overhang < 0 {
		for i := 0; i < dSeq.Overhang*-1; i++ {
			fmt.Print(" ")
		}
	}
	fmt.Println(dSeq.Crick)
}

func NewFromWatsonStrand(watson string, geometry constants.SequenceGeometry) *Dseq {
	return &Dseq{
		Watson:   watson,
		Crick:    transform.Complement(watson),
		Overhang: 0,
		Geometry: geometry,
	}
}

func NewDseq(watson, crick string, overhang int, geometry constants.SequenceGeometry) *Dseq {
	return &Dseq{
		Watson:   watson,
		Crick:    crick,
		Overhang: overhang,
		Geometry: geometry,
	}
}

type Cutter interface {
	GetNextRecognitionSite(sequence string, offset int, isCircular bool) []enzyme.RecognitionSiteResult
}

/*
Cut the Dseq with the provided enzyme.

Takes either a single enzyme or a batch of enzymes (RestrictionBatch).

Returns a slice of Dseqs that represent the fragments of the Dseq after cutting with the enzyme.

If the enzyme is circular, the sequence is treated as circular and the
function will return the next recognition site, even if it spans the
beginning and end of the sequence.
*/
func (dSeq *Dseq) Cut(enzyme Cutter) []Dseq {
	fragments := make([]Dseq, 0)

	nextSearchStart := 0
	lastOverhang := 0

	lastCrickCutIndex := 0
	lastWatsonCutIndex := 0

	for {
		results := enzyme.GetNextRecognitionSite(
			dSeq.Watson,
			nextSearchStart,
			dSeq.Geometry == constants.Circular,
		)
		if results == nil {
			break
		}

		var watsonCutIndex, crickCutIndex int

		if len(results) > 1 {
			fmt.Println("Warning: multiple recognition sites found. Only the first site will be used.")
		}
		result := results[0]
		enzyme := result.Enzyme

		if result.Strand == constants.Watson {
			watsonCutIndex = result.RecognitionSiteIndex + enzyme.FivePrimeCutSite
			crickCutIndex = result.RecognitionSiteIndex + enzyme.ThreePrimeCutSite

			if watsonCutIndex >= len(dSeq.Watson) {
				watsonCutIndex = len(dSeq.Watson)
			}

			if crickCutIndex >= len(dSeq.Crick) {
				crickCutIndex = len(dSeq.Crick)
			}

			fragment := Dseq{
				Watson:   dSeq.Watson[lastWatsonCutIndex:watsonCutIndex],
				Crick:    dSeq.Crick[lastCrickCutIndex:crickCutIndex],
				Overhang: lastOverhang,
				Geometry: constants.Linear,
			}
			fragments = append(fragments, fragment)

			lastOverhang = (crickCutIndex - watsonCutIndex) * -1
			nextSearchStart = watsonCutIndex

		} else {
			watsonCutIndex = result.RecognitionSiteIndex + enzyme.Length - enzyme.ThreePrimeCutSite
			crickCutIndex = result.RecognitionSiteIndex + enzyme.Length - enzyme.FivePrimeCutSite
			overhang := watsonCutIndex - crickCutIndex

			fragment := Dseq{
				Watson:   dSeq.Watson[lastWatsonCutIndex:watsonCutIndex],
				Crick:    dSeq.Crick[lastCrickCutIndex:crickCutIndex],
				Overhang: overhang,
				Geometry: constants.Linear,
			}
			fragments = append(fragments, fragment)

			nextSearchStart = result.RecognitionSiteIndex + 1

		}

		lastCrickCutIndex = crickCutIndex
		lastWatsonCutIndex = watsonCutIndex
	}

	// Add the last fragment
	fragment := Dseq{
		Watson:   dSeq.Watson[lastWatsonCutIndex:],
		Crick:    dSeq.Crick[lastCrickCutIndex:],
		Overhang: lastOverhang,
		Geometry: dSeq.Geometry,
	}
	fragments = append(fragments, fragment)

	return fragments
}
