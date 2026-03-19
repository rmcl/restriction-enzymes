package enzyme

import (
	"sort"
	"strings"

	"github.com/rmcl/restriction-enzymes/constants"
)

type RestrictionBatch struct {
	Enzymes            []Enzyme
	fwdSiteToEnzyme    map[string][]*Enzyme
	revSiteToEnzyme    map[string][]*Enzyme
	orderedEnzymes     []*Enzyme
	exactFwdByLength   map[int]map[string][]*Enzyme
	exactRevByLength   map[int]map[string][]*Enzyme
	exactSiteLengths   []int
	maxIndexedSiteLen  int
	nonExactEnzymes    []*Enzyme
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
	restrictionBatch.resetIndex()
}

// Add all the enzymes from another restriction batch to this one.
func (restrictionBatch *RestrictionBatch) AddBatch(batch RestrictionBatch) {
	restrictionBatch.Enzymes = append(restrictionBatch.Enzymes, batch.Enzymes...)
	restrictionBatch.resetIndex()
}

// Return a mapping of enzyme name to a list of sites in the sequence it cuts.
// This returns the position of the cut site on the watson strand to mirror
// Biopython's interface.
//
// If isCircular is true, the sequence is treated as circular and the
// function will return a recognition site, even if it spans the
// beginning and end of the sequence.
func (enzymeBatch *RestrictionBatch) Search(sequence string, isCircular bool) (map[string][]int, error) {
	enzymeBatch.buildIndex()

	watsonCutSitesByEnzyme := make(map[string][]int)
	normalizedSequence := strings.ToLower(sequence)
	exactSearchSequence := buildSearchSequence(normalizedSequence, enzymeBatch.maxIndexedSiteLen, isCircular)

	for position := 0; position < len(normalizedSequence); position++ {
		results := enzymeBatch.getExactResultsForPosition(position, exactSearchSequence)
		for _, result := range results {
			enzyme := result.Enzyme
			if enzymeCutSites, ok := watsonCutSitesByEnzyme[enzyme.Name]; ok {
				enzymeCutSites = append(enzymeCutSites, result.WatsonCutIndex)
				watsonCutSitesByEnzyme[enzyme.Name] = enzymeCutSites
			} else {
				watsonCutSitesByEnzyme[enzyme.Name] = []int{result.WatsonCutIndex}
			}
		}
	}

	for _, restrictionEnzyme := range enzymeBatch.nonExactEnzymes {
		positions, err := restrictionEnzyme.Search(normalizedSequence, isCircular)
		if err != nil {
			return nil, err
		}
		if len(positions) == 0 {
			continue
		}
		watsonCutSitesByEnzyme[restrictionEnzyme.Name] = append(
			watsonCutSitesByEnzyme[restrictionEnzyme.Name],
			positions...,
		)
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

	restrictionBatch.buildIndex()

	sequence = strings.ToLower(sequence)
	exactSearchSequence := buildSearchSequence(sequence, restrictionBatch.maxIndexedSiteLen, isCircular)

	for position := offset; position < len(sequence); position++ {
		results := restrictionBatch.getExactResultsForPosition(position, exactSearchSequence)
		results = append(results, restrictionBatch.getNonExactResultsForPosition(position, sequence, isCircular)...)
		if len(results) > 0 {
			return results
		}
	}

	return nil
}

func (restrictionBatch *RestrictionBatch) getExactResultsForPosition(
	matchStart int,
	sequence string,
) []RecognitionSiteResult {
	results := []RecognitionSiteResult{}
	matchedEnzymeNames := map[string]bool{}
	for _, siteLength := range restrictionBatch.exactSiteLengths {
		if matchStart+siteLength > len(sequence) {
			continue
		}

		siteSequence := sequence[matchStart : matchStart+siteLength]

		for _, candidate := range restrictionBatch.exactFwdByLength[siteLength][siteSequence] {
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

		for _, candidate := range restrictionBatch.exactRevByLength[siteLength][siteSequence] {
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

func (restrictionBatch *RestrictionBatch) getNonExactResultsForPosition(
	matchStart int,
	sequence string,
	isCircular bool,
) []RecognitionSiteResult {
	results := []RecognitionSiteResult{}
	matchedEnzymeNames := map[string]bool{}

	for _, candidate := range restrictionBatch.nonExactEnzymes {
		remainingSequence := buildSearchSequence(sequence[matchStart:], len(candidate.Site), isCircular)

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

func (restrictionBatch *RestrictionBatch) resetIndex() {
	restrictionBatch.fwdSiteToEnzyme = nil
	restrictionBatch.revSiteToEnzyme = nil
	restrictionBatch.orderedEnzymes = nil
	restrictionBatch.exactFwdByLength = nil
	restrictionBatch.exactRevByLength = nil
	restrictionBatch.exactSiteLengths = nil
	restrictionBatch.maxIndexedSiteLen = 0
	restrictionBatch.nonExactEnzymes = nil
}

func (restrictionBatch *RestrictionBatch) buildIndex() {
	if restrictionBatch.orderedEnzymes != nil {
		return
	}

	restrictionBatch.fwdSiteToEnzyme = make(map[string][]*Enzyme)
	restrictionBatch.revSiteToEnzyme = make(map[string][]*Enzyme)
	restrictionBatch.orderedEnzymes = []*Enzyme{}
	restrictionBatch.exactFwdByLength = make(map[int]map[string][]*Enzyme)
	restrictionBatch.exactRevByLength = make(map[int]map[string][]*Enzyme)
	restrictionBatch.exactSiteLengths = []int{}
	restrictionBatch.maxIndexedSiteLen = 0
	restrictionBatch.nonExactEnzymes = []*Enzyme{}

	if len(restrictionBatch.Enzymes) == 0 {
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

	for i := range enzymes {
		enzyme := &enzymes[i]
		restrictionBatch.orderedEnzymes = append(restrictionBatch.orderedEnzymes, enzyme)

		// Unfortunately, we need to strip the first 4 characters from the
		// the regular expression pattern because it is "(?i)"
		fwdSite := strings.ToLower(enzyme.RegexpFor.String()[4:])
		revSite := strings.ToLower(enzyme.RegexpRev.String()[4:])

		if restrictionBatch.fwdSiteToEnzyme[fwdSite] == nil {
			restrictionBatch.fwdSiteToEnzyme[fwdSite] = []*Enzyme{}
		}
		if restrictionBatch.revSiteToEnzyme[revSite] == nil {
			restrictionBatch.revSiteToEnzyme[revSite] = []*Enzyme{}
		}

		restrictionBatch.fwdSiteToEnzyme[fwdSite] = append(restrictionBatch.fwdSiteToEnzyme[fwdSite], enzyme)
		restrictionBatch.revSiteToEnzyme[revSite] = append(restrictionBatch.revSiteToEnzyme[revSite], enzyme)

		if isExactRecognitionSite(fwdSite) && isExactRecognitionSite(revSite) {
			if restrictionBatch.exactFwdByLength[enzyme.Length] == nil {
				restrictionBatch.exactFwdByLength[enzyme.Length] = map[string][]*Enzyme{}
				restrictionBatch.exactRevByLength[enzyme.Length] = map[string][]*Enzyme{}
				restrictionBatch.exactSiteLengths = append(restrictionBatch.exactSiteLengths, enzyme.Length)
			}
			restrictionBatch.exactFwdByLength[enzyme.Length][fwdSite] = append(restrictionBatch.exactFwdByLength[enzyme.Length][fwdSite], enzyme)
			restrictionBatch.exactRevByLength[enzyme.Length][revSite] = append(restrictionBatch.exactRevByLength[enzyme.Length][revSite], enzyme)
			if enzyme.Length > restrictionBatch.maxIndexedSiteLen {
				restrictionBatch.maxIndexedSiteLen = enzyme.Length
			}
		} else {
			restrictionBatch.nonExactEnzymes = append(restrictionBatch.nonExactEnzymes, enzyme)
		}
	}

	sort.Slice(restrictionBatch.exactSiteLengths, func(i, j int) bool {
		return restrictionBatch.exactSiteLengths[i] > restrictionBatch.exactSiteLengths[j]
	})
}

func buildSearchSequence(sequence string, maxSiteLen int, isCircular bool) string {
	if !isCircular || maxSiteLen <= 1 || len(sequence) == 0 {
		return sequence
	}

	requiredExtension := maxSiteLen - 1
	if requiredExtension <= 0 {
		return sequence
	}

	extendedSequence := sequence
	for len(extendedSequence) < len(sequence)+requiredExtension {
		extendedSequence += sequence
	}

	return extendedSequence[:len(sequence)+requiredExtension]
}

func isExactRecognitionSite(site string) bool {
	for _, nucleotide := range site {
		if nucleotide != 'a' && nucleotide != 'c' && nucleotide != 'g' && nucleotide != 't' {
			return false
		}
	}
	return true
}
