package enzyme

import (
	"regexp"
	"strings"
	"testing"

	"github.com/rmcl/restriction-enzymes/constants"
)

func TestGetNextRecognitionSite(t *testing.T) {
	seq1 := "NNGGTCTCNCACANNNNCTCTNGAGACCNN"

	enzyme := FIXTURES["BsaI"]

	result := enzyme.GetNextRecognitionSite(seq1, 0, true)
	if result[0].RecognitionSiteIndex != 2 {
		t.Errorf("Expected 2, got %d", result[0].RecognitionSiteIndex)
		return
	}
	if result[0].Strand != constants.Watson {
		t.Errorf("Expected Watson, got %s", result[0].Strand)
		return
	}

	result = enzyme.GetNextRecognitionSite(seq1, 3, true)
	if result[0].RecognitionSiteIndex != 22 {
		t.Errorf("Expected 22, got %d", result[0].RecognitionSiteIndex)
		return
	}
	if result[0].Strand != constants.Crick {
		t.Errorf("Expected Watson, got %s", result[0].Strand)
		return
	}
}

func TestGetNextRecognitionSiteWithLowerCase(t *testing.T) {
	seq1 := "NNGGTCTCNCACANNNNCTCTNGAGACCNN"
	seq1 = strings.ToLower(seq1)

	enzyme := FIXTURES["BsaI"]

	results := enzyme.GetNextRecognitionSite(seq1, 0, true)
	if results[0].RecognitionSiteIndex != 2 {
		t.Errorf("Expected 2, got %d", results[0].RecognitionSiteIndex)
		return
	}
	if results[0].Strand != constants.Watson {
		t.Errorf("Expected Watson, got %s", results[0].Strand)
		return
	}

	results = enzyme.GetNextRecognitionSite(seq1, 3, true)
	if results[0].RecognitionSiteIndex != 22 {
		t.Errorf("Expected 22, got %d", results[0].RecognitionSiteIndex)
		return
	}
	if results[0].Strand != constants.Crick {
		t.Errorf("Expected Watson, got %s", results[0].Strand)
		return
	}
}

func TestCutSeqLongSeqNewEnzyme(t *testing.T) {

	enzyme := FIXTURES["EcoRI"]

	results := enzyme.GetNextRecognitionSite(EXAMPLE_SEQUENCE_1, 0, true)

	if results[0].RecognitionSiteIndex != 3807 {
		t.Errorf("Expected 3807, got %d", results[0].RecognitionSiteIndex)
		return
	}
	if results[0].Strand != constants.Watson {
		t.Errorf("Expected crick, got %s", results[0].Strand)
		return
	}
}

func TestCutSeqWithWrappedSite(t *testing.T) {
	seq1 := "TCNNNCACANNNNCTCTNGAGNNACCNNGGTC"

	enzyme := FIXTURES["BsaI"]

	results := enzyme.GetNextRecognitionSite(seq1, 0, true)
	if results[0].RecognitionSiteIndex != 28 {
		t.Errorf("Expected 28, got %d", results[0].RecognitionSiteIndex)
		return
	}
	if results[0].Strand != constants.Watson {
		t.Errorf("Expected Watson, got %s", results[0].Strand)
		return
	}
}

func TestCutSeqRandomLongSequence(t *testing.T) {

	enzyme := FIXTURES["BamHI"]

	results := enzyme.GetNextRecognitionSite(EXAMPLE_SEQUENCE_2, 0, true)
	if results[0].RecognitionSiteIndex != 824 {
		t.Errorf("Expected 824, got %d", results[0].RecognitionSiteIndex)
		return
	}
	if results[0].Strand != constants.Watson {
		t.Errorf("Expected Watson, got %s", results[0].Strand)
		return
	}
}

func TestBsaIExampleGetRecognitionSite(t *testing.T) {
	enzyme := FIXTURES["BsaI"]

	results := enzyme.GetNextRecognitionSite(EXAMPLE_SEQUENCE_3, 0, true)
	if results[0].RecognitionSiteIndex != 165 {
		t.Errorf("Expected 165, got %d", results[0].RecognitionSiteIndex)
		return
	}
	if results[0].Strand != constants.Crick {
		t.Errorf("Expected Crick, got %s", results[0].Strand)
		return
	}
}

func TestBsaISearch(t *testing.T) {
	enzyme := FIXTURES["BsaI"]

	results, err := enzyme.Search(EXAMPLE_SEQUENCE_3, false)
	if err != nil {
		t.Errorf("Error: %s", err)
		return
	}

	if len(results) != 2 {
		t.Errorf("Expected 2, got %d", len(results))
		return
	}

	if results[0] != 160 {
		t.Errorf("Expected 160, got %d", results[0])
		return
	}

	if results[1] != 590 {
		t.Errorf("Expected 590, got %d", results[1])
		return
	}
}

var FIXTURES = map[string]Enzyme{
	"BsaI":  {Name: "BsaI", Site: "GGTCTC", Length: 6, Substrate: "DNA", RegexpFor: regexp.MustCompile("(?i)GGTCTC"), RegexpRev: regexp.MustCompile("(?i)GAGACC"), OverhangLength: 0, OverhangSequence: "", NumberOfCuts: 2, CutType: "sticky", FivePrimeCutSite: 7, ThreePrimeCutSite: 11, FivePrimeCutSite2: 0, ThreePrimeCutSite2: 0, RebaseId: 313, InactivationTemperature: 0, OptimalTemperature: 0, Uri: "https://identifiers.org/rebase:313", References: []string{"Flodman K., Xu S.-Y.;", "Fomenkov A.;", "Kong H., Chen Z.;", "Morgan R.D.;", "Xu S.-Y.;", "Zhu Z., Xu S.-Y.;", "Zhu Z., Xu S.-Y.;"}},
	"EcoRI": {Name: "EcoRI", Site: "GAATTC", Length: 6, Substrate: "DNA", RegexpFor: regexp.MustCompile("(?i)GAATTC"), RegexpRev: regexp.MustCompile("(?i)GAATTC"), OverhangLength: 0, OverhangSequence: "", NumberOfCuts: 2, CutType: "sticky", FivePrimeCutSite: 1, ThreePrimeCutSite: 5, FivePrimeCutSite2: 0, ThreePrimeCutSite2: 0, RebaseId: 993, InactivationTemperature: 0, OptimalTemperature: 0, Uri: "https://identifiers.org/rebase:993", References: []string{"Albertsen H.M., Le Paslier D., Abderrahim H., Dausset J., Cann H., ", "Cohen D.;", "Dugaiczyk A., Hedgpeth J., Boyer H.W., Goodman H.M.;", "Flodman K., Xu S.-Y.;", "Forrow S., Lee M., Souhami R.L., Hartley J.A.;", "Greene P.J., Betlach M.C., Boyer H.W., Goodman H.M.;", "Hedgpeth J., Goodman H.M., Boyer H.W.;", "Newman A.K., Rubin R.A., Kim S.H., Modrich P.;", "Pingoud A., Alves J., Fliess A., Geiger R., Rueter T., Wolfes H.;", "Tanaka M.;", "Winkler F.K.;"}},
	"BamHI": {Name: "BamHI", Site: "GGATCC", Length: 6, Substrate: "DNA", RegexpFor: regexp.MustCompile("(?i)GGATCC"), RegexpRev: regexp.MustCompile("(?i)GGATCC"), OverhangLength: 0, OverhangSequence: "", NumberOfCuts: 2, CutType: "sticky", FivePrimeCutSite: 1, ThreePrimeCutSite: 5, FivePrimeCutSite2: 0, ThreePrimeCutSite2: 0, RebaseId: 185, InactivationTemperature: 0, OptimalTemperature: 0, Uri: "https://identifiers.org/rebase:185", References: []string{"Brooks J.E., Nathan P.D., Landry D., Sznyter L.A., Waite-Rees P., ", "Ives C.L., Moran L.S., Slatko B.E., Benner J.S.;", "Endo M., Majima T.;", "Flodman K., Xu S.-Y.;", "Fomenkov A.;", "Fomenkov A., Anton B.P., Vincze T., Roberts R.J.;", "Majima T., Endo M.;", "Roberts R.J., Wilson G.A., Young F.E.;", "Usami S., Kurimura H., Kino K., Kamigaki K., Kirimura K.;", "Wilson G.A., Young F.E.;"}},
}
var EXAMPLE_SEQUENCE_1 = "ATGACATAACCGTATTACCGCCATGCATTAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCATAGCCCATATATGGAGTTCCGCGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCAACGACCCCCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCCAATAGGGACTTTCCATTGACGTCAATGGGTGGAGTATTTACGGTAAACTGCCCACTTGGCAGTACATCAAGTGTATCATATGCCAAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCCTGGCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCGGTTTGACTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCTGGTTTAGTGAACCGTCAGATCCGCTAGTCGACGGTACCTCGGCGATGGCTTTTCCGCCGCGGCGACGGCTGCGCCTCGGTCCCCGCGGCCTCCCGCTTCTTCTCTCGGGACTCCTGCTACCTCTGTGCCGCGCCTTCAACCTAGACGTGGACAGTCCTGCCGAGTACTCTGGCCCCGAGGGAAGTTACTTCGGCTTCGCCGTGGATTTCTTCGTGCCCAGCGCGTCTTCCCGGATGTTTCTTCTCGTGGGAGCTCCCAAAGCAAACACCACCCAGCCTGGGATTGTGGAAGGAGGGCAGGTCCTCAAATGTGACTGGTCTTCTACCCGCCGGTGCCAGCCAATTGAATTTGATGCAACAGGCAATAGAGATTATGCCAAGGATGATCCATTGGAATTTAAGTCCCATCAGTGGTTTGGAGCATCTGTGAGGTCGAAACAGGATAAAATTTTGGCCTGTGCCCCATTGTACCATTGGAGAACTGAGATGAAACAGGAGCGAGAGCCTGTTGGAACATGCTTTCTTCAAGATGGAACAAAGACTGTTGAGTATGCTCCATGTAGATCACAAGATATTGATGCTGATGGACAGGGATTTTGTCAAGGAGGATTCAGCATTGATTTTACTAAAGCTGACAGAGTACTTCTTGGTGGTCCTGGTAGCTTTTATTGGCAAGGTCAGCTTATTTCGGATCAAGTGGCAGAAATCGTATCTAAATACGACCCCAATGTTTACAGCATCAAGTATAATAACCAATTAGCAACTCGGACTGCACAAGCTATTTTTGATGACAGCTATTTGGGTTATTCTGTGGCTGTCGGAGATTTCAATGGTGATGGCATAGATGACTTTGTTTCAGGAGTTCCAAGAGCAGCAAGGACTTTGGGAATGGTTTATATTTATGATGGGAAGAACATGTCCTCCTTATACAATTTTACTGGCGAGCAGATGGCTGCATATTTCGGATTTTCTGTAGCTGCCACTGACATTAATGGAGATGATTATGCAGATGTGTTTATTGGAGCACCTCTCTTCATGGATCGTGGCTCTGATGGCAAACTCCAAGAGGTGGGGCAGGTCTCAGTGTCTCTACAGAGAGCTTCAGGAGACTTCCAGACGACAAAGCTGAATGGATTTGAGGTCTTTGCACGGTTTGGCAGTGCCATAGCTCCTTTGGGAGATCTGGACCAGGATGGTTTCAATGATATTGCAATTGCTGCTCCATATGGGGGTGAAGATAAAAAAGGAATTGTTTATATCTTCAATGGAAGATCAACAGGCTTGAACGCAGTCCCATCTCAAATCCTTGAAGGGCAGTGGGCTGCTCGAAGCATGCCACCAAGCTTTGGCTATTCAATGAAAGGAGCCACAGATATAGACAAAAATGGATATCCAGACTTAATTGTAGGAGCTTTTGGTGTAGATCGAGCTATCTTATACAGGGCCAGACCAGTTATCACTGTAAATGCTGGTCTTGAAGTGTACCCTAGCATTTTAAATCAAGACAATAAAACCTGCTCACTGCCTGGAACAGCTCTCAAAGTTTCCTGTTTTAATGTTAGGTTCTGCTTAAAGGCAGATGGCAAAGGAGTACTTCCCAGGAAACTTAATTTCCAGGTGGAACTTCTTTTGGATAAACTCAAGCAAAAGGGAGCAATTCGACGAGCACTGTTTCTCTACAGCAGGTCCCCAAGTCACTCCAAGAACATGACTATTTCAAGGGGGGGACTGATGCAGTGTGAGGAATTGATAGCGTATCTGCGGGATGAATCTGAATTTAGAGACAAACTCACTCCAATTACTATTTTTATGGAATATCGGTTGGATTATAGAACAGCTGCTGATACAACAGGCTTGCAACCCATTCTTAACCAGTTCACGCCTGCTAACATTAGTCGACAGGCTCACATTCTACTTGACTGTGGTGAAGACAATGTCTGTAAACCCAAGCTGGAAGTTTCTGTAGATAGTGATCAAAAGAAGATCTATATTGGGGATGACAACCCTCTGACATTGATTGTTAAGGCTCAGAATCAAGGAGAAGGTGCCTACGAAGCTGAGCTCATCGTTTCCATTCCACTGCAGGCTGATTTCATCGGGGTTGTCCGAAACAATGAAGCCTTAGCAAGACTTTCCTGTGCATTTAAGACAGAAAACCAAACTCGCCAGGTGGTATGTGACCTTGGAAACCCAATGAAGGCTGGAACTCAACTCTTAGCTGGTCTTCGTTTCAGTGTGCACCAGCAGTCAGAGATGGATACTTCTGTGAAATTTGACTTACAAATCCAAAGCTCAAATCTATTTGACAAAGTAAGCCCAGTTGTATCTCACAAAGTTGATCTTGCTGTTTTAGCTGCAGTTGAGATAAGAGGAGTCTCGAGTCCTGATCATATCTTTCTTCCGATTCCAAACTGGGAGCACAAGGAGAACCCTGAGACTGAAGAAGATGTTGGGCCAGTTGTTCAGCACATCTATGAGCTGAGAAACAATGGTCCAAGTTCATTCAGCAAGGCAATGCTCCATCTTCAGTGGCCTTACAAATATAATAATAACACTCTGTTGTATATCCTTCATTATGATATTGATGGACCAATGAACTGCACTTCAGATATGGAGATCAACCCTTTGAGAATTAAGATCTCATCTTTGCAAACAACTGAAAAGAATGACACGGTTGCCGGGCAAGGTGAGCGGGACCATCTCATCACTAAGCGGGATCTTGCCCTCAGTGAAGGAGATATTCACACTTTGGGTTGTGGAGTTGCTCAGTGCTTGAAGATTGTCTGCCAAGTTGGGAGATTAGACAGAGGAAAGAGTGCAATCTTGTACGTAAAGTCATTACTGTGGACTGAGACTTTTATGAATAAAGAAAATCAGAATCATTCCTATTCTCTGAAGTCGTCTGCTTCATTTAATGTCATAGAGTTTCCTTATAAGAATCTTCCAATTGAGGATATCACCAACTCCACATTGGTTACCACTAATGTCACCTGGGGCATTCAGCCAGCGCCCATGCCTGTGCCTGTGTGGGTGATCATTTTAGCAGTTCTAGCAGGATTGTTGCTACTGGCTGTTTTGGTATTTGTAATGTACAGGATGGGCTTTTTTAAACGGGTCCGGCCACCTCAAGAAGAACAAGAAAGGGAGCAGCTTCAACCTCATGAAAATGGTGAAGGAAACTCAGAAACTCCGGGATCTCGAGCTCAAGCTTCGAATTCTGCAGTCGACGGTACCGCGGGCCCGGGATCCCCACCGGTCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCTTGACCTACGGCGTGCAGTGCTTCGCCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAAGGTCTATATCACCGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGACCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCAAGCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCGGCCGCGACTCTAGATCATAATCAGCCATACCACATTTGTAGAGGTTTTACTTGCTTTAAAAAACCTCCCACACCTCCCCCTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTAAGGCGTAAATTGTAAGCGTTAATATTTTGTTAAAATTCGCGTTAAATTTTTGTTAAATCAGCTCATTTTTTAACCAATAGGCCGAAATCGGCAAAATCCCTTATAAATCAAAAGAATAGACCGAGATAGGGTTGAGTGTTGTTCCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGGCCCACTACGTGAACCATCACCCTAATCAAGTTTTTTGGGGTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAGAAAGCGAAAGGAGCGGGCGCTAGGGCGCTGGCAAGTGTAGCGGTCACGCTGCGCGTAACCACCACACCCGCCGCGCTTAATGCGCCGCTACAGGGCGCGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTCCTGAGGCGGAAAGAACCAGCTGTGGAATGTGTGTCAGTTAGGGTGTGGAAAGTCCCCAGGCTCCCCAGCAGGCAGAAGTATGCAAAGCATGCATCTCAATTAGTCAGCAACCAGGTGTGGAAAGTCCCCAGGCTCCCCAGCAGGCAGAAGTATGCAAAGCATGCATCTCAATTAGTCAGCAACCATAGTCCCGCCCCTAACTCCGCCCATCCCGCCCCTAACTCCGCCCAGTTCCGCCCATTCTCCGCCCCATGGCTGACTAATTTTTTTTATTTATGCAGAGGCCGAGGCCGCCTCGGCCTCTGAGCTATTCCAGAAGTAGTGAGGAGGCTTTTTTGGAGGCCTAGGCTTTTGCAAAGATCGATCAAGAGACAGGATGAGGATCGTTTCGCATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCAGCGCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTGCAAGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTGCTCGACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAGGATCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATGCGGCGGCTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGCATCGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAAGAGCATCAGGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGAGCATGCCCGACGGCGAGGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAATGGCCGCTTTTCTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGACATAGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTCCTCGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTTGACGAGTTCTTCTGAGCGGGACTCTGGGGTTCGAAATGACCGACCAAGCGACGCCCAACCTGCCATCACGAGATTTCGATTCCACCGCCGCCTTCTATGAAAGGTTGGGCTTCGGAATCGTTTTCCGGGACGCCGGCTGGATGATCCTCCAGCGCGGGGATCTCATGCTGGAGTTCTTCGCCCACCCTAGGGGGAGGCTAACTGAAACACGGAAGGAGACAATACCGGAAGGAACCCGCGCTATGACGGCAATAAAAAGACAGAATAAAACGCACGGTGTTGGGTCGTTTGTTCATAAACGCGGGGTTCGGTCCCAGGGCTGGCACTCTGTCGATACCCCACCGAGACCCCATTGGGGCCAATACGCCCGCGTTTCTTCCTTTTCCCCACCCCACCCCCCAAGTTCGGGTGAAGGCCCAGGGCTCGCAGCCAACGTCGGGGCGGCAGGCCCTGCCATAGCCTCAGGTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGAA"
var EXAMPLE_SEQUENCE_2 = "ttacggggtcattagttcatagcccatatatggagttccgcgttacataacttacggtaaatggcccgcctggctgaccgcccaacgacccccgcccattgacgtcaataatgacgtatgttcccatagtaacgccaatagggactttccattgacgtcaatgggtggagtatttacggtaaactgcccacttggcagtacatcaagtgtatcatatgccaagtacgccccctattgacgtcaatgacggtaaatggcccgcctggcattatgcccagtacatgaccttatgggactttcctacttggcagtacatctacgtattagtcatcgctattaccatggtgatgcggttttggcagtacatcaatgggcgtggatagcggtttgactcacggggatttccaagtctccaccccattgacgtcaatgggagtttgttttggcaccaaaatcaacgggactttccaaaatgtcgtaacaactccgccccattgacgcaaatgggcggtaggcgtgtacggtgggaggtctatataagcagagctggtttagtgaaccgtcagatccgctagcatgaggcttcgggagccgctcctgagcggcagcgccgcgatgccaggcgcgtccctacagcgggcctgccgcctgctcgtggccgtctgcgctctgcaccttggcgtcaccctcgtttactacctggctggccgcgacctgagccgcctgccccaactggtcggagtctccacaccgctgcagggcggctcgaacagtgccgccgccatcgggcagtcctccggggagctccggaccggaggggccaaggatccaccggtcgccaccatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttgacctacggcgtgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaaggtctatatcaccgccgacaagcagaagaacggcatcaaggtgaacttcaagacccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccaagctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaagcggccgcgactctagatcataatcagccataccacatttgtagaggttttacttgctttaaaaaacctcccacacctccccctgaacctgaaacataaaatgaatgcaattgttgttgttaacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatcttaaggcgtaaattgtaagcgttaatattttgttaaaattcgcgttaaatttttgttaaatcagctcattttttaaccaataggccgaaatcggcaaaatcccttataaatcaaaagaatagaccgagatagggttgagtgttgttccagtttggaacaagagtccactattaaagaacgtggactccaacgtcaaagggcgaaaaaccgtctatcagggcgatggcccactacgtgaaccatcaccctaatcaagttttttggggtcgaggtgccgtaaagcactaaatcggaaccctaaagggagcccccgatttagagcttgacggggaaagccggcgaacgtggcgagaaaggaagggaagaaagcgaaaggagcgggcgctagggcgctggcaagtgtagcggtcacgctgcgcgtaaccaccacacccgccgcgcttaatgcgccgctacagggcgcgtcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtcctgaggcggaaagaaccagctgtggaatgtgtgtcagttagggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccaggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccatagtcccgcccctaactccgcccatcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcctaggcttttgcaaagatcgatcaagagacaggatgaggatcgtttcgcatgattgaacaagatggattgcacgcaggttctccggccgcttgggtggagaggctattcggctatgactgggcacaacagacaatcggctgctctgatgccgccgtgttccggctgtcagcgcaggggcgcccggttctttttgtcaagaccgacctgtccggtgccctgaatgaactgcaagacgaggcagcgcggctatcgtggctggccacgacgggcgttccttgcgcagctgtgctcgacgttgtcactgaagcgggaagggactggctgctattgggcgaagtgccggggcaggatctcctgtcatctcaccttgctcctgccgagaaagtatccatcatggctgatgcaatgcggcggctgcatacgcttgatccggctacctgcccattcgaccaccaagcgaaacatcgcatcgagcgagcacgtactcggatggaagccggtcttgtcgatcaggatgatctggacgaagagcatcaggggctcgcgccagccgaactgttcgccaggctcaaggcgagcatgcccgacggcgaggatctcgtcgtgacccatggcgatgcctgcttgccgaatatcatggtggaaaatggccgcttttctggattcatcgactgtggccggctgggtgtggcggaccgctatcaggacatagcgttggctacccgtgatattgctgaagagcttggcggcgaatgggctgaccgcttcctcgtgctttacggtatcgccgctcccgattcgcagcgcatcgccttctatcgccttcttgacgagttcttctgagcgggactctggggttcgaaatgaccgaccaagcgacgcccaacctgccatcacgagatttcgattccaccgccgccttctatgaaaggttgggcttcggaatcgttttccgggacgccggctggatgatcctccagcgcggggatctcatgctggagttcttcgcccaccctagggggaggctaactgaaacacggaaggagacaataccggaaggaacccgcgctatgacggcaataaaaagacagaataaaacgcacggtgttgggtcgtttgttcataaacgcggggttcggtcccagggctggcactctgtcgataccccaccgagaccccattggggccaatacgcccgcgtttcttccttttccccaccccaccccccaagttcgggtgaaggcccagggctcgcagccaacgtcggggcggcaggccctgccatagcctcaggttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgccatgcattagttattaatagtaatcaa"
var EXAMPLE_SEQUENCE_3 = "ccggcgtagaggatcgagatcgatctcgatcccgcgaaattaatacgactcactataggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacatatgagccatcatcatcaccatcatggctgctcaggaggagaccgaaagcgtattacagtgacagttgacagcgacagctatcagttgctcaaggcatatatgatgtcaatatctccggtctggtatgcacaaccaagaatgaagcccatgcagtttaaggtttacacctataaaagagagagccgttatcgtctgtttgtggatgtacagagtgatattattgacacgcccgggcgacggatggtgatccccctggccagtgcacgtctgctgtcagataaagtctcccgtgaactttacccggtggtgcatatcggggatgaaagctggcgcatgatgaccaccgatatggccagtgtgccggtttccgttatcggggaagaagtggctgatctcagccaccgcgaaaatgacatcaaaaacgccattaacctgatgttctggggaatatgaacggtctcgttcctaatgagatccggctgctaacaaagcccgaaaggaagctgagttggctgctgccaccgctgagcaataactagcataaccccttggggcctctaaacgggtcttgaggggttttttgctgaaaggaggaactatatccggattggcgaatgggacgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgtttacaatttcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgaattaattcttagaaaaactcatcgagcatcaaatgaaactgcaatttattcatatcaggattatcaataccatatttttgaaaaagccgtttctgtaatgaaggagaaaactcaccgaggcagttccataggatggcaagatcctggtatcggtctgcgattccgactcgtccaacatcaatacaacctattaatttcccctcgtcaaaaataaggttatcaagtgagaaatcaccatgagtgacgactgaatccggtgagaatggcaaaagtttatgcatttctttccagacttgttcaacaggccagccattacgctcgtcatcaaaatcactcgcatcaaccaaaccgttattcattcgtgattgcgcctgagcgagacgaaatacgcgatcgctgttaaaaggacaattacaaacaggaatcgaatgcaaccggcgcaggaacactgccagcgcatcaacaatattttcacctgaatcaggatattcttctaatacctggaatgctgttttcccggggatcgcagtggtgagtaaccatgcatcatcaggagtacggataaaatgcttgatggtcggaagaggcataaattccgtcagccagtttagtctgaccatctcatctgtaacatcattggcaacgctacctttgccatgtttcagaaacaactctggcgcatcgggcttcccatacaatcgatagattgtcgcacctgattgcccgacattatcgcgagcccatttatacccatataaatcagcatccatgttggaatttaatcgcggcctagagcaagacgtttcccgttgaatatggctcataacaccccttgtattactgtttatgtaagcagacagttttattgttcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagtatacactccgctatcgctacgtgactgggtcatggctgcgccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgaggcagctgcggtaaagctcatcagcgtggtcgtgaagcgattcacagatgtctgcctgttcatccgcgtccagctcgttgagtttctccagaagcgttaatgtctggcttctgataaagcgggccatgttaagggcggttttttcctgtttggtcactgatgcctccgtgtaagggggatttctgttcatgggggtaatgataccgatgaaacgagagaggatgctcacgatacgggttactgatgatgaacatgcccggttactggaacgttgtgagggtaaacaactggcggtatggatgcggcgggaccagagaaaaatcactcagggtcaatgccagcgcttcgttaatacagatgtaggtgttccacagggtagccagcagcatcctgcgatgcagatccggaacataatggtgcagggcgctgacttccgcgtttccagactttacgaaacacggaaaccgaagaccattcatgttgttgctcaggtcgcagacgttttgcagcagcagtcgcttcacgttcgctcgcgtatcggtgattcattctgctaaccagtaaggcaaccccgccagcctagccgggtcctcaacgacaggagcacgatcatgcgcacccgtggggccgccatgccggcgataatggcctgcttctcgccgaaacgtttggtggcgggaccagtgacgaaggcttgagcgagggcgtgcaagattccgaataccgcaagcgacaggccgatcatcgtcgcgctccagcgaaagcggtcctcgccgaaaatgacccagagcgctgccggcacctgtcctacgagttgcatgataaagaagacagtcataagtgcggcgacgatagtcatgccccgcgcccaccggaaggagctgactgggttgaaggctctcaagggcatcggtcgagatcccggtgcctaatgagtgagctaacttacattaattgcgttgcgctcactgcccgctttccagtcgggaaacctgtcgtgccagctgcattaatgaatcggccaacgcgcggggagaggcggtttgcgtattgggcgccagggtggtttttcttttcaccagtgagacgggcaacagctgattgcccttcaccgcctggccctgagagagttgcagcaagcggtccacgctggtttgccccagcaggcgaaaatcctgtttgatggtggttaacggcgggatataacatgagctgtcttcggtatcgtcgtatcccactaccgagatgtccgcaccaacgcgcagcccggactcggtaatggcgcgcattgcgcccagcgccatctgatcgttggcaaccagcatcgcagtgggaacgatgccctcattcagcatttgcatggtttgttgaaaaccggacatggcactccagtcgccttcccgttccgctatcggctgaatttgattgcgagtgagatatttatgccagccagccagacgcagacgcgccgagacagaacttaatgggcccgctaacagcgcgatttgctggtgacccaatgcgaccagatgctccacgcccagtcgcgtaccgtcttcatgggagaaaataatactgttgatgggtgtctggtcagagacatcaagaaataacgccggaacattagtgcaggcagcttccacagcaatggcatcctggtcatccagcggatagttaatgatcagcccactgacgcgttgcgcgagaagattgtgcaccgccgctttacaggcttcgacgccgcttcgttctaccatcgacaccaccacgctggcacccagttgatcggcgcgagatttaatcgccgcgacaatttgcgacggcgcgtgcagggccagactggaggtggcaacgccaatcagcaacgactgtttgcccgccagttgttgtgccacgcggttgggaatgtaattcagctccgccatcgccgcttccactttttcccgcgttttcgcagaaacgtggctggcctggttcaccacgcgggaaacggtctgataagagacaccggcatactctgcgacatcgtataacgttactggtttcacattcaccaccctgaattgactctcttccgggcgctatcatgccataccgcgaaaggttttgcgccattcgatggtgtccgggatctcgacgctctcccttatgcgactcctgcattaggaagcagcccagtagtaggttgaggccgttgagcaccgccgccgcaaggaatggtgcatgcaaggagatggcgcccaacagtcccccggccacggggcctgccaccatacccacgccgaaacaagcgctcatgagcccgaagtggcgagcccgatcttccccatcggtgatgtcggcgatataggcgccagcaaccgcacctgtggcgccggtgatgccggccacgatgcgt"
