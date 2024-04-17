package enzyme

import (
	"testing"

	"github.com/rmcl/restriction-enzymes/constants"
)

func TestGetNextRecognitionSiteForBatch(unittest *testing.T) {

	batch := NewRestrictionBatch(
		FIXTURES["BsaI"],
		FIXTURES["EcoRI"],
	)
	sequence := "AATAGACAGAATTCGATTCACCAGAGGTCTCATAGACAAACCAGAGAAAAAAAA"

	position, enzyme, strand := batch.GetNextRecognitionSite(sequence, 0, false)
	if (position != 8) || enzyme.Name != "EcoRI" || (strand != constants.Watson) {
		unittest.Errorf("Expected 8, nil, EcoRI, got %d, %v, %s", position, enzyme, strand)
	}

	position, enzyme, strand = batch.GetNextRecognitionSite(sequence, 9, false)
	if (position != 25) || enzyme.Name != "BsaI" || (strand != constants.Watson) {
		unittest.Errorf("Expected 25, nil, BsaI, got %d, %v, %s", position, enzyme, strand)
	}
}

func TestAddEnzymeToBatch(unittest *testing.T) {

	batch := NewRestrictionBatch(
		FIXTURES["BsaI"],
		FIXTURES["EcoRI"],
	)

	if len(batch.Enzymes) != 2 {
		unittest.Errorf("Expected 2 enzymes, got %d", len(batch.Enzymes))
	}

	batch.Add(FIXTURES["BamHI"])

	if len(batch.Enzymes) != 3 {
		unittest.Errorf("Expected 3 enzymes, got %d", len(batch.Enzymes))
	}

}
