package sequence

import (
	"fmt"

	"github.com/bebop/poly/transform"
	"github.com/rmcl/restriction-enzymes/enzyme"
)

type Strand string

const (
	Watson Strand = "watson"
	Crick  Strand = "crick"
)

type SequenceGeometry string

const (
	Linear   SequenceGeometry = "linear"
	Circular SequenceGeometry = "circular"
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
	Geometry SequenceGeometry
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

func NewFromWatsonStrand(watson string, geometry SequenceGeometry) *Dseq {
	return &Dseq{
		Watson:   watson,
		Crick:    transform.Complement(watson),
		Overhang: 0,
		Geometry: geometry,
	}
}

func NewDseq(watson, crick string, overhang int, geometry SequenceGeometry) *Dseq {
	return &Dseq{
		Watson:   watson,
		Crick:    crick,
		Overhang: overhang,
		Geometry: geometry,
	}
}

func GetNextRecognitionSite(
	sequence string,
	offset int,
	enzyme enzyme.Enzyme,
	isCircular bool,
) (int, Strand) {

	sequence = sequence[offset:]

	// If the enzyme is circular, we need to support the case where the
	// recognition site spans the beginning and end of the sequence.
	// To do this, we append the sequence to itself of length site - 1.
	if isCircular {
		sequence += sequence[:len(enzyme.Site)-1]
	}

	var watsonMatchIndex int
	var crickMatchIndex int

	watsonMatch := enzyme.RegexpFor.FindStringIndex(sequence)
	if watsonMatch != nil {
		watsonMatchIndex = watsonMatch[0] + offset
	} else {
		watsonMatchIndex = -1
	}

	crickMatch := enzyme.RegexpRev.FindStringIndex(sequence)
	if crickMatch != nil {
		crickMatchIndex = crickMatch[0] + offset
	} else {
		crickMatchIndex = -1
	}

	if watsonMatchIndex > -1 && crickMatchIndex > -1 {

		if watsonMatchIndex < crickMatchIndex {
			return watsonMatchIndex, Watson
		} else {
			return crickMatchIndex, Crick
		}

	} else if watsonMatchIndex > -1 {
		return watsonMatchIndex, Watson
	} else {
		return crickMatchIndex, Crick
	}

}

func (dSeq *Dseq) Cut(enzyme enzyme.Enzyme) []Dseq {
	fragments := make([]Dseq, 0)

	nextSearchStart := 0
	lastSiteIndex := 0
	lastOverhang := 0

	for {
		siteIndex, strand := GetNextRecognitionSite(dSeq.Watson, nextSearchStart, enzyme, dSeq.Geometry == Circular)
		if siteIndex == -1 {
			break
		}

		if strand == Watson {
			watsonCutIndex := siteIndex + enzyme.FivePrimeCutSite
			crickCutIndex := siteIndex + enzyme.ThreePrimeCutSite

			fragment := Dseq{
				Watson:   dSeq.Watson[lastSiteIndex:watsonCutIndex],
				Crick:    dSeq.Crick[lastSiteIndex:crickCutIndex],
				Overhang: lastOverhang,
				Geometry: dSeq.Geometry,
			}
			fragments = append(fragments, fragment)

			lastOverhang = (crickCutIndex - watsonCutIndex) * -1
			lastSiteIndex = watsonCutIndex
			nextSearchStart = lastSiteIndex
			fmt.Println("Watson", lastOverhang)

		} else {
			fmt.Println("Crick", lastOverhang, siteIndex, enzyme.Length, enzyme.FivePrimeCutSite, enzyme.ThreePrimeCutSite)

			watsonCutIndex := siteIndex + enzyme.Length - enzyme.ThreePrimeCutSite
			crickCutIndex := siteIndex + enzyme.Length - enzyme.FivePrimeCutSite
			overhang := watsonCutIndex - crickCutIndex

			fragment := Dseq{
				Watson:   dSeq.Watson[lastSiteIndex:watsonCutIndex],
				Crick:    dSeq.Crick[lastSiteIndex-overhang : crickCutIndex],
				Overhang: overhang,
				Geometry: dSeq.Geometry,
			}
			fragments = append(fragments, fragment)

			nextSearchStart = siteIndex + 1
			lastSiteIndex = watsonCutIndex

		}
	}

	// Add the last fragment
	fragment := Dseq{
		Watson:   dSeq.Watson[lastSiteIndex:],
		Crick:    dSeq.Crick[lastSiteIndex-lastOverhang:],
		Overhang: lastOverhang,
		Geometry: dSeq.Geometry,
	}
	fragments = append(fragments, fragment)

	return fragments
}
