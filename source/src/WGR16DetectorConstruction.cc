//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: WGR16DetectorConstruction.cc 77656 2013-11-27 08:52:57Z gcosmo $
//
/// \file WGR16DetectorConstruction.cc
/// \brief Implementation of the WGR16DetectorConstruction class
#include "WGR16PMTSD.hh"
#include "WGR16DetectorConstruction.hh"
#include "WGR16MagneticField.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
// #include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4GenericTrap.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UserLimits.hh"
#include "G4PVParameterised.hh"
#include "G4ThreeVector.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"
#include "G4VisExtent.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "geomdefs.hh"

#include <cmath>
#include <stdio.h>
#include <float.h>

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal WGR16MagneticField* WGR16DetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* WGR16DetectorConstruction::fFieldMgr = 0;
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16DetectorConstruction::WGR16DetectorConstruction()
: G4VUserDetectorConstruction(), fMessenger(0), fVisAttributes(),fScoringVolume(0)
{
	DefineCommands();
}


WGR16DetectorConstruction::~WGR16DetectorConstruction()
{
	delete fMessenger;
	
	for (G4int i=0; i<G4int(fVisAttributes.size()); ++i)
	{
		delete fVisAttributes[i];
	}
}


G4VPhysicalVolume* WGR16DetectorConstruction::Construct()
{
	G4bool checkOverlaps = false;

	//////////////////////
	// -- meterials -- //
	/////////////////////
	ConstructMaterials();
	G4Material* vac = G4Material::GetMaterial("G4_Galactic");
	G4Material* vac_PMTHouse = G4Material::GetMaterial("G4_Galactic");
	G4Material* cu  = G4Material::GetMaterial("G4_Cu");
	// const G4double cu_radlen = cu->GetRadlen(); // -- radiation length -- //

	G4double z, a, density, fractionmass;
	G4int ncomponents, natoms;
	G4String symbol;
	G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
	G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
	G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
	G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
	G4Element* F  = new G4Element("Fluorine",symbol="F" , z= 9., a= 18.9984*g/mole);
	// G4Element* Si = new G4Element("Silicon" ,symbol="Si", z= 14., a= 28.09*g/mole);/

	// -- for PMT Cathod -- //	
	G4Material* Al 
	= new G4Material("Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);

	// G4double PMTT = 1*mm;

	// -- Photocathod property -- //
	G4MaterialPropertiesTable* mpPMTPC = this->MaterialPropertyTable_PMTPC();
	Al->SetMaterialPropertiesTable( mpPMTPC );
	G4OpticalSurface* OpSurf_PMTPC = new G4OpticalSurface("OpSurf_PMTPC",glisur,polished,dielectric_metal);
	OpSurf_PMTPC->SetMaterialPropertiesTable(mpPMTPC);

	G4MaterialPropertiesTable* mpPMTHouse = this->MaterialPropertyTable_PMTHouse();
	vac_PMTHouse->SetMaterialPropertiesTable( mpPMTHouse );
	G4OpticalSurface* OpSurf_PMTHouse = new G4OpticalSurface("OpSurf_PMTHouse",unified,polished,dielectric_metal);
	OpSurf_PMTHouse->SetMaterialPropertiesTable(mpPMTHouse);




	// new G4Material("Lead"     , z=82., a=207.19*g/mole, density=11.35*g/cm3);
	// new G4Material("Copper"   , z=29., a=63.546*g/mole, density=8.96*g/cm3);

	// -- for PMT Glass -- //
	G4Material* Glass = new G4Material("Glass", density=1.032*g/cm3,2);
	Glass->AddElement(C,91.533*perCent);
	Glass->AddElement(H,8.467*perCent);

	G4MaterialPropertiesTable* mpGlass = this->MaterialPropertyTable_Glass();
	Glass->SetMaterialPropertiesTable(mpGlass);


	// -- for scintillation fiber core -- //
	G4Material* polystyrene
	= new G4Material("Polystyrene",density= 1.05*g/cm3, ncomponents=2);
	polystyrene->AddElement(C, natoms=8);
	polystyrene->AddElement(H, natoms=8);


	// -- for cladding (scintillation fibers) -- //
	G4Material* pmma_clad
	= new G4Material("PMMA_Clad",density= 1.19*g/cm3, ncomponents=3);
	pmma_clad->AddElement(C, natoms=5);
	pmma_clad->AddElement(H, natoms=8);
	pmma_clad->AddElement(O, natoms=2);


	// -- for Cerenkov fiber core -- //
	G4Material* pmma 
	= new G4Material("PMMA",density= 1.19*g/cm3, ncomponents=3);
	pmma->AddElement(C, natoms=5);
	pmma->AddElement(H, natoms=8);
	pmma->AddElement(O, natoms=2);

	G4MaterialPropertiesTable* mpPMMA = this->MaterialPropertyTable_PMMA();
	pmma->SetMaterialPropertiesTable(mpPMMA);


	// -- for cladding (Cerenkov fibers) -- //
	G4Material* fluorinatedPolymer 
	= new G4Material("Fluorinated_Polymer", density= 1.43*g/cm3, ncomponents=2);
	fluorinatedPolymer->AddElement(C,2);
	fluorinatedPolymer->AddElement(F,2);

	G4MaterialPropertiesTable* mpFS = this->MaterialPropertyTable_FS();
	fluorinatedPolymer->SetMaterialPropertiesTable(mpFS);


	// -- Air -- //
	G4Material* Air 
	= new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
	Air->AddElement(N, fractionmass=0.7);
	Air->AddElement(O, fractionmass=0.3);

	G4MaterialPropertiesTable* mpAir = this->MaterialPropertyTable_Air();
	Air->SetMaterialPropertiesTable(mpAir);

	// G4MaterialPropertiesTable* mpPS
	// G4MaterialPropertiesTable* mpPMTPC;

	/////////////////
	// -- world -- //
	/////////////////
	G4VSolid* worldSolid 
	= new G4Box("worldBox",10.*m,10.*m,10.*m);
	G4LogicalVolume* worldLogical
	= new G4LogicalVolume(worldSolid,vac,"worldLogical");
	G4VPhysicalVolume* worldPhysical
	= new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0, false,0,checkOverlaps);

	//////////////////
	// -- Cu box -- //
	//////////////////
	G4double radius = 1.8*m; // -- size of inner tracker: it should be empty space to avoid any overlap with tracker system -- //
	G4double nTower_PhiDir = 283;
//	 G4double nTower_PhiDir = 30;

	G4double pi = 3.14159265358979323846;
	G4double dPhi = (2*pi) / nTower_PhiDir;
	G4double half_dPhi = 0.5*dPhi;
	G4double CuLen_PhiDir = 2*radius*std::tan(half_dPhi);
	G4double CuLen_EtaDir = 39.9648624417329649818*mm;
	G4double CuLen_H = 2.5*m;
	cout<<CuLen_PhiDir<<" "<< CuLen_EtaDir<<endl;

	//finding set of etas and saving it Eta
	G4double Eta_Max=atan(2.5*m/radius);
	G4double Eta[52];
	Eta[0]=0;
	const G4int nTower_EtaDir=51;
	for (G4int i=1;i<=nTower_EtaDir;++i){
		Eta[i]=this->WGR16DetectorConstruction::Solve_Eq(Eta[i-1],Eta[i-1],Eta_Max,radius,CuLen_EtaDir);
	}


	bool DrawOneUnitTower = false;
	if( DrawOneUnitTower )
	{
		if( nTower_PhiDir != 283 ) cout << "nTower_PhiDir = " << nTower_PhiDir << " should be 283 for correct geometry" << endl;
		nTower_PhiDir = 1;
		CuLen_EtaDir = 40*mm;
		CuLen_PhiDir = CuLen_EtaDir;
		CuLen_H = 2.5*m;
	}

	cout << "[Cu] (PhiDir, EtaDir, Height) = (" << CuLen_PhiDir << ", " << CuLen_EtaDir << ", " << CuLen_H << ")" << endl;

	G4Box* CuBox
	= new G4Box("CuBox", CuLen_EtaDir/2.0, CuLen_PhiDir/2.0, CuLen_H/2.0);
	G4LogicalVolume *CuLogical
	= new G4LogicalVolume(CuBox, cu, "CuLogical");




	//////////////////
	// -- Fibers -- //
	//////////////////
	// -- Materials for Cerenkov fiber -- //
	G4Material *clad_C_Material = fluorinatedPolymer;
	G4Material *core_C_Material = pmma;

	// -- Materials for Scintillation fiber -- //
	G4Material *clad_S_Material = pmma_clad;
	G4Material *core_S_Material = polystyrene;

	// // -- Material for PMT glass -- //
	// G4Material *Glass_Material = Glass;

	// // -- Material for PMT Photocathod -- //
	// G4Material *PMTPC_Material = Al;

	// -- fibre parameters -- //
	G4double clad_C_rMin = 0*mm;
	G4double clad_C_rMax = 0.50*mm;
	G4double clad_C_Dz   = 2.5*m;
	G4double clad_C_Sphi = 0.;
	G4double clad_C_Dphi = 2.*M_PI;

	G4double core_C_rMin = 0.*mm;
	G4double core_C_rMax = 0.49*mm;
	G4double core_C_Dz   = 2.5*m;
	G4double core_C_Sphi = 0.;
	G4double core_C_Dphi = 2.*M_PI;

	// G4double clad_S_rMin = 0.485*mm;
	// G4double clad_S_rMax = 0.50*mm;
	// G4double clad_S_Dz   = 2.5*m;
	// G4double clad_S_Sphi = 0.;
	// G4double clad_S_Dphi = 2.*M_PI;

	G4double core_S_rMin = 0.*mm;
	G4double core_S_rMax = 0.485*mm;
	G4double core_S_Dz   = 2.5*m;
	G4double core_S_Sphi = 0.;
	G4double core_S_Dphi = 2.*M_PI;

	G4double dist_btwCore = 1.5*mm;

	bool Do_test = false;
	if( nTower_PhiDir == 10 ) Do_test = true;
	if( Do_test )
	{
		if( nTower_PhiDir != 10 )
			cout << "This test setting may not work for this # tower = " << nTower_PhiDir << endl;
		dist_btwCore = 250*mm; // -- when # tower = 10 -- //
		clad_C_rMax = dist_btwCore / 3.0;
		core_C_rMax = clad_C_rMax * 0.98;
		core_S_rMax = clad_C_rMax * 0.98;
	}

	const G4int nFiber_PhiDir = floor( CuLen_PhiDir / dist_btwCore ) - 1;
	const G4int nFiber_EtaDir = floor( CuLen_EtaDir / dist_btwCore ) - 1;
	G4double dist_edge_PhiDir = ( CuLen_PhiDir - (nFiber_PhiDir-1)*dist_btwCore ) / 2.0;
	G4double dist_edge_EtaDir = ( CuLen_EtaDir - (nFiber_EtaDir-1)*dist_btwCore ) / 2.0;

	cout << "nFiber_PhiDir: " << nFiber_PhiDir << endl;
	cout << "nFiber_EtaDir: " << nFiber_EtaDir << endl;
	cout << "dist_edge_PhiDir: " << dist_edge_PhiDir << ", dist_edge_EtaDir: " << dist_edge_EtaDir << endl;

	// -- Solid -- //
	G4VSolid* fiberClad = new G4Tubs("fiberClad", 0.*mm, clad_C_rMax, CuLen_H/2., clad_C_Sphi, clad_C_Dphi); // -- S is the same -- //
	G4VSolid* fiberCoreC = new G4Tubs("fiberCoreC", 0.*mm, core_C_rMax, CuLen_H/2., core_C_Sphi, core_C_Dphi);
	G4VSolid* fiberCoreS = new G4Tubs("fiberCoreS", 0.*mm, core_S_rMax, CuLen_H/2., core_S_Sphi, core_S_Dphi);

	// -- PMT House, glass and PhotoCathod -- //
	G4double PMTHouseLen_H = 3*mm;
	G4double PMTHouseLen_EtaDir = CuLen_EtaDir;
	G4double PMTHouseLen_PhiDir = CuLen_PhiDir;

	G4double PMTGlassLen_H = 2*mm;
	G4double PMTGlassLen_EtaDir	= PMTHouseLen_EtaDir;
	G4double PMTGlassLen_PhiDir = PMTHouseLen_PhiDir;

	G4double PMTPCLen_H = 1*mm;
	G4double PMTPCLen_EtaDir = PMTHouseLen_EtaDir;
	G4double PMTPCLen_PhiDir = PMTHouseLen_PhiDir;

	if( Do_test )
	{
		PMTHouseLen_H *= 10;
		PMTGlassLen_H *= 10;
		PMTPCLen_H *= 10;
	}

	G4Box* PMTHouseBox 
	= new G4Box("PMTHouseBox", PMTHouseLen_EtaDir/2.0, PMTHouseLen_PhiDir/2., PMTHouseLen_H/2.0);
	G4LogicalVolume *PMTHouseBox_Logic
	= new G4LogicalVolume(PMTHouseBox, vac_PMTHouse, "PMTHouseBox_Logic");

	G4Box* PMTGlassBox 
	= new G4Box("PMTGlassBox", PMTGlassLen_EtaDir/2.0, PMTGlassLen_PhiDir/2., PMTGlassLen_H/2.0);
	G4LogicalVolume *PMTGlassBox_Logic
	= new G4LogicalVolume(PMTGlassBox, Glass, "PMTGlassBox_Logic");

	G4Box* PMTPCBox 
	= new G4Box("PMTPCBox", PMTPCLen_EtaDir/2.0, PMTPCLen_PhiDir/2., PMTPCLen_H/2.0);
	this->PMTPCBox_Logic
	= new G4LogicalVolume(PMTPCBox, Al, "PMTPCBox_Logic");

	
	// new G4LogicalSkinSurface("SkinSurf_PMTHouse", PMTHouseBox_Logic, OpSurf_PMTHouse);
	new G4LogicalSkinSurface("SkinSurf_PMTPC", PMTPCBox_Logic, OpSurf_PMTPC);
	G4int fiberNumber = 0;
	// -- iteration for eta direction -- //
	for(G4int i_barrel=-51; i_barrel<=51; i_barrel++)
	{
		G4double dEta;
		G4double dEta_pre;
		if (i_barrel<0) {
			dEta=-Eta[-i_barrel];
			dEta_pre=-Eta[-i_barrel-1];
		}
		else {
			dEta=Eta[i_barrel];
			if (i_barrel!=0) dEta_pre=Eta[i_barrel-1];
		}

		
		//////////////////////////////////////////////////////////////////////////////////
		// -- Logical volume of Cu trapezoid for phi direction (triangle) -- //
		////////////////////////////////////////////////////////////////////////////////
		G4double Angle_btwTower_cos = pow(cos(dEta),2)*cos(dPhi) + pow(sin(dEta),2);
		G4double Angle_btwTower_tan = sqrt(1-pow(Angle_btwTower_cos,2))/Angle_btwTower_cos;

		G4double CuTrdAngle_face = acos(pow(sin(dEta),2)*cos(dPhi) + pow(cos(dEta),2))/2.0;
		G4double diffOfPhiDir = 2*CuLen_EtaDir * sin(CuTrdAngle_face);
		G4double CuTrdLen_PhiDirShort = CuLen_H*Angle_btwTower_tan;
		G4double CuTrdLen_PhiDirLong = CuTrdLen_PhiDirShort + diffOfPhiDir;
		
		vector<G4TwoVector> Trap_Vertices;

		Trap_Vertices.push_back(G4TwoVector(-CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
		Trap_Vertices.push_back(G4TwoVector(-CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
		Trap_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
		Trap_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0 - diffOfPhiDir*sin(CuTrdAngle_face),diffOfPhiDir*cos(CuTrdAngle_face) + CuLen_PhiDir/2.0));
		Trap_Vertices.push_back(G4TwoVector(-CuLen_EtaDir/2.0 - CuTrdLen_PhiDirShort*sin(CuTrdAngle_face),CuTrdLen_PhiDirShort*cos(CuTrdAngle_face) + CuLen_PhiDir/2.0));
		Trap_Vertices.push_back(G4TwoVector(-CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
		Trap_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
		Trap_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0 - CuTrdLen_PhiDirLong*sin(CuTrdAngle_face),CuTrdLen_PhiDirLong*cos(CuTrdAngle_face) + CuLen_PhiDir/2.0));

		G4GenericTrap *CuTrap
		= new G4GenericTrap("CuTrap", CuLen_H/2.0,Trap_Vertices);
		G4LogicalVolume* CuTrapLogical
		= new G4LogicalVolume(CuTrap, cu, "CuTrapLogical");

		//////////////////////////////////////////////////////////////////////////////////////////////
		// -- Logical volume of Cu trapezoid for Eta direction for rectangular part -- //
		/////////////////////////////////////////////////////////////////////////////////////////////
		G4double TrapAngle;
		G4double CuTrapLen_EtaDir;
		G4GenericTrap *CuTrapEta;
		G4LogicalVolume* CuTrapEta_Logical;
		if (i_barrel!=0){
			TrapAngle = dEta-dEta_pre;
			CuTrapLen_EtaDir = CuLen_H*fabs(tan(TrapAngle));
			vector<G4TwoVector> TrapEta_Vertices;
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,-CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,-CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,-CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuTrapLen_EtaDir + CuLen_EtaDir/2.0,-CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuTrapLen_EtaDir + CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
			CuTrapEta
			= new G4GenericTrap("CuTrapEta", CuLen_H/2.0,TrapEta_Vertices);
			CuTrapEta_Logical
			= new G4LogicalVolume(CuTrapEta, cu, "CuTrapEtaLogical");
			
		}

		//////////////////////////////////////////////////////////////////////////////////////////////
		// -- Logical volume of Cu trapezoid for Eta direction for triangular part -- //
		/////////////////////////////////////////////////////////////////////////////////////////////
		G4double diffOfPhiDir_Trap;
		G4double CuTrdLen_PhiDirShort_Pre;
		G4GenericTrap *CuTrapTrdEta;
		G4LogicalVolume* CuTrapTrdEta_Logical;
		if (i_barrel!=0){
			diffOfPhiDir_Trap = 2*CuTrapLen_EtaDir*sin(CuTrdAngle_face);
			G4double CuTrdLen_PhiDirShort_Pre = CuTrdLen_PhiDirLong + diffOfPhiDir_Trap;
			vector<G4TwoVector> TrapEta_Vertices;
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0 - diffOfPhiDir*sin(CuTrdAngle_face),diffOfPhiDir*cos(CuTrdAngle_face) + CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0 - diffOfPhiDir*sin(CuTrdAngle_face),diffOfPhiDir*cos(CuTrdAngle_face) + CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0 - CuTrdLen_PhiDirLong*sin(CuTrdAngle_face),CuTrdLen_PhiDirLong*cos(CuTrdAngle_face) + CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuLen_EtaDir/2.0,CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuTrapLen_EtaDir + CuLen_EtaDir/2.0, CuLen_PhiDir/2.0));
			TrapEta_Vertices.push_back(G4TwoVector(CuTrapLen_EtaDir - CuTrdLen_PhiDirShort_Pre*sin(CuTrdAngle_face) + CuLen_EtaDir/2.0,CuTrdLen_PhiDirShort_Pre*cos(CuTrdAngle_face) + CuLen_PhiDir/2.0));
			CuTrapTrdEta
			= new G4GenericTrap("CuTrapTrdEta", CuLen_H/2.0,TrapEta_Vertices);
			CuTrapTrdEta_Logical
			= new G4LogicalVolume(CuTrapTrdEta, cu, "CuTrapEtaLogical");
		}
		// -- iteration for phi direction -- //
		for(G4int i_cu=0; i_cu<2; i_cu++)
		{

			////////////////////
			// -- Cu box -- //
			////////////////////
			G4double phi = i_cu*dPhi;
			G4RotationMatrix rotM  = G4RotationMatrix();
			if (i_barrel<0) rotM.rotateZ(180*deg);
			rotM.rotateY(90*degree-dEta);
			rotM.rotateZ(phi);
			G4double Phi_move = radius + 0.5*CuLen_H*cos(dEta) + 0.5*CuLen_EtaDir*fabs(sin(dEta));
			G4double Eta_move = radius*tan(dEta) + 0.5*CuLen_H*sin(dEta) + 0.5*CuLen_EtaDir*fabs(sin(dEta))*tan(dEta);
			G4ThreeVector Trans_Vector = G4ThreeVector(Phi_move*cos(phi),Phi_move*sin(phi),Eta_move);
			G4Transform3D transform = G4Transform3D(rotM,Trans_Vector);

			new G4PVPlacement(transform, CuLogical, to_string(i_barrel)+" "+to_string(i_cu) + " CuPhysical", worldLogical, false, i_cu, checkOverlaps );
			////////////////////////////
			// -- Cu trapezoid -- //
			//////////////////////////
			new G4PVPlacement(transform, CuTrapLogical, to_string(i_barrel)+" "+to_string(i_cu)+" CuTrapPhysical", worldLogical, false, i_cu, checkOverlaps );
			/////////////////////////////////////////////////////
			//-- Cu trapezoidEta for rectangular part & triangular part -- //
			//////////////////////////////////////////////////////
			if (i_barrel!=0){
				new G4PVPlacement(transform, CuTrapEta_Logical, to_string(i_barrel)+" CuTrapEtaPhysical", worldLogical, false, i_cu, checkOverlaps );
				new G4PVPlacement(transform, CuTrapTrdEta_Logical, to_string(i_barrel)+" CuTrapTrdEtaPhysical", worldLogical, false, i_cu, checkOverlaps );
			}
			
			// -- PMTs -- //
			// G4ThreeVector position_PMTHouse = (radius + CuLen_H + 0.5*PMTHouseLen_H)*Unit_Z;
			// G4Transform3D transform_PMTHouse = G4Transform3D(rotM,position_PMTHouse);
			// new G4PVPlacement(transform_PMTHouse, PMTHouseBox_Logic, "PMTHouseBox_Phys", worldLogical, false, i_cu, checkOverlaps );

			// new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*PMTHouseLen_H + 0.5*PMTGlassLen_H), PMTGlassBox_Logic, "PMTGlassBox_Phys", PMTHouseBox_Logic, false, i_cu, checkOverlaps );
			// new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*PMTHouseLen_H - 0.5*PMTPCLen_H), this->PMTPCBox_Logic, "PMTPCBox_Phys", PMTHouseBox_Logic, false, i_cu, checkOverlaps );

//			////////////////////////////
//			// -- Cu box fiber -- //
//			///////////////////////////
//			 G4int i_total = 0;
//			 G4int fiberC_total = 0;
//			 G4int fiberS_total = 0;
//			 for(G4int i_EtaDir=0; i_EtaDir<nFiber_EtaDir; i_EtaDir++)
//			 {
//			 	G4double x_EtaDir = ((-1)*CuLen_EtaDir / 2.0) + dist_edge_EtaDir + dist_btwCore*i_EtaDir;
//			 	for(G4int i_PhiDir=0; i_PhiDir<nFiber_PhiDir; i_PhiDir++)
//			 	{
//			 		i_total++;
//			 		G4double x_PhiDir = ((-1)*CuLen_PhiDir / 2.0) + dist_edge_PhiDir + dist_btwCore*i_PhiDir;
//
//			 		// -- cladding: same shape for both C and S fiber -- //
//			 		G4VSolid* FiberClad_ith
//			 		= new G4IntersectionSolid("fiberClad", CuBox, fiberClad, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//
//			 		G4LogicalVolume *FiberClad_Logic_ith
//			 		= new G4LogicalVolume(FiberClad_ith, clad_C_Material, "FiberClad_Logic");
//
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberClad_Logic_ith, "FiberClad_Phys", CuLogical, false, i_total, checkOverlaps);
//
//			 		// -- Cores -- //
//			 		G4VSolid* FiberCore_ith;
//			 		G4LogicalVolume *FiberCore_Logic_ith;
//			 		// -- s, c, s, c, ... -- //
//			 		// -- c, s, c, s, ... -- //
//			 		bool isFiberC = this->IsFiberC(i_EtaDir, i_PhiDir);
//			 		if( isFiberC )
//			 		{
//			 			++fiberC_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuBox, fiberCoreC, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_C_Material, "FiberCore_Logic");
//			 		}
//			 		else
//			 		{
//			 			++fiberS_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuBox, fiberCoreS, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_S_Material, "FiberCore_Logic");
//			 		}
//
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberCore_Logic_ith, "FiberCore_Phys", FiberClad_Logic_ith, false, i_total, checkOverlaps);
//
//			 		G4VisAttributes* visAttr = new G4VisAttributes();
//			 		if( isFiberC )
//			 			visAttr->SetColour( G4Colour(0.0,0.0,1.0) ); // -- blue -- //
//			 		else
//			 			visAttr->SetColour( G4Colour(1.0,1.0,0.0) );  // -- yellow -- //
//			 		visAttr->SetForceSolid(true);
//			 		visAttr->SetVisibility(true);
//			 		FiberCore_Logic_ith->SetVisAttributes(visAttr);
//			 	}
//			 }
//			 cout<<i_barrel<<" CuBox"<<endl;
//			 cout<<"Cerenkov fiber #: "<<fiberC_total<<"		Scintillator fiber #: "<<fiberS_total<<"		total fiber #: "<<i_total<<endl;

//			/////////////////////////////////
//			// -- Cu trapezoid fiber-- //
//			////////////////////////////////
//			
//			// -- tri_fibers -- //
//			G4double CuTrdLen_EtaDir = CuLen_EtaDir*cos(CuTrdAngle_face);
//			const G4int tri_nFiber_EtaDir = floor( CuTrdLen_EtaDir / dist_btwCore ) - 1;
//			G4double tri_dist_edge_EtaDir = ( CuTrdLen_EtaDir - (tri_nFiber_EtaDir-1)*dist_btwCore ) / 2.0;
//			const G4int tri_nFiber_PhiDirOuter = floor( CuTrdLen_PhiDirLong / dist_btwCore ) - 1;
//			G4double tri_dist_edge_PhiDirOuter = ( CuTrdLen_PhiDirLong - (tri_nFiber_PhiDirOuter-1)*dist_btwCore ) / 2.0;
//
//			G4int tri_i_total = 0;
//			G4int tri_fiberC_total = 0;
//			G4int tri_fiberS_total = 0;
//			for(G4int i_EtaDir=0; i_EtaDir< tri_nFiber_EtaDir; i_EtaDir++)
//			{
//				
//				G4double CuTrdLen_PhiDir_Fiber
//				= CuTrdLen_PhiDirShort + 2*(tri_dist_edge_EtaDir + dist_btwCore*i_EtaDir)*tan(CuTrdAngle_face);
//				const G4int tri_nFiber_PhiDir = floor( CuTrdLen_PhiDir_Fiber / dist_btwCore ) - 1;
//
//				 for(G4int i_PhiDir=0; i_PhiDir<tri_nFiber_PhiDir; i_PhiDir++)
//				 {
//				 	tri_i_total++;
//				 	G4double x_EtaDir
//				 	= ((-1)*CuLen_EtaDir / 2.0) + (tri_dist_edge_EtaDir + dist_btwCore*i_EtaDir)/cos(CuTrdAngle_face) - (tri_dist_edge_PhiDirOuter + dist_btwCore*i_PhiDir)*sin(CuTrdAngle_face);
//				 	G4double x_PhiDir
//					= CuLen_PhiDir/2.0 + (tri_dist_edge_PhiDirOuter + dist_btwCore*i_PhiDir)*cos(CuTrdAngle_face);
//				 	
//				 	G4ThreeVector trans_vector=G4ThreeVector(x_EtaDir, x_PhiDir, 0);
//				 	
//				 	// -- cladding: same shape for both C and S fiber -- //
//			 		G4VSolid *FiberClad_ith
//					= new G4IntersectionSolid("fiberClad", CuTrap, fiberClad, 0, trans_vector);
//			 		
//			 		G4LogicalVolume *FiberClad_Logic_ith
//			 		= new G4LogicalVolume(FiberClad_ith, clad_C_Material, "FiberClad_Logic");
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberClad_Logic_ith, "FiberClad_Phys", CuTrapLogical, false, tri_i_total, checkOverlaps);
//				 	
//			 		// -- Cores -- //
//			 		G4VSolid* FiberCore_ith;
//			 		G4LogicalVolume *FiberCore_Logic_ith;
//				 	// -- s, c, s, c, ... -- //
//			 		// -- c, s, c, s, ... -- //
//			 		bool isFiberC = this->IsFiberC(i_EtaDir, i_PhiDir);
//			 		if( isFiberC )
//			 		{
//			 			++tri_fiberC_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuTrap, fiberCoreC, 0, trans_vector);
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_C_Material, "FiberCore_Logic");
//			 		}
//			 		else
//			 		{
//			 			++tri_fiberS_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuTrap, fiberCoreS, 0, trans_vector);
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_S_Material, "FiberCore_Logic");
//			 		}
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberCore_Logic_ith, "FiberCore_Phys", FiberClad_Logic_ith, false, tri_i_total, checkOverlaps);
//				 		G4VisAttributes* visAttr1 = new G4VisAttributes();
//			 		if( isFiberC ){
//			 			visAttr1->SetColour( G4Colour(0.0,0.0,1.0) ); // -- blue -- //
//			 		}
//			 		else{
//			 			visAttr1->SetColour( G4Colour(1.0,1.0,0.0) );  // -- yellow -- //
//			 		}
//			 		visAttr1->SetForceSolid(true);
//			 		visAttr1->SetVisibility(true);
//			 		FiberCore_Logic_ith->SetVisAttributes(visAttr1);
//			 	}
//			 }
//		cout<<i_barrel<<" CuTrapPhi"<<endl;
//		cout<<"Cerenkov fiber #: "<<tri_fiberC_total<<"		Scintillator fiber #: "<<tri_fiberS_total<<"		total fiber #: "<<tri_i_total<<endl;
//
//			
//			///////////////////////////////////////////////////////////
//			//-- Cu trapezoidEta for rectangular part fiber-- //
//			///////////////////////////////////////////////////////////
//
//			// -- Trdfibers for rectangular part-- //
//
//			const G4int Trd_nFiber_EtaDir = floor( CuTrapLen_EtaDir / dist_btwCore) - 1;
//			G4double Trd_dist_edge_EtaDir = ( CuTrapLen_EtaDir - (Trd_nFiber_EtaDir-1)*dist_btwCore)/ 2.0;
//			G4int TrapEta_i_total = 0;
//			G4int TrapEta_fiberC_total = 0;
//			G4int TrapEta_fiberS_total = 0;
//
//			for(G4int i_EtaDir=0; i_EtaDir<Trd_nFiber_EtaDir; i_EtaDir++)
//			{
//				G4double x_EtaDir = CuLen_PhiDir / 2.0 + Trd_dist_edge_EtaDir + dist_btwCore*i_EtaDir;
//			 	for(G4int i_PhiDir=0; i_PhiDir<nFiber_PhiDir; i_PhiDir++)
//			 	{
//			 		TrapEta_i_total++;
//			 		G4double x_PhiDir = ((-1)*CuLen_PhiDir / 2.0) + dist_edge_PhiDir + dist_btwCore*i_PhiDir;
//
//			 		// -- cladding: same shape for both C and S fiber -- //
//			 		G4VSolid* FiberClad_ith
//			 		= new G4IntersectionSolid("fiberClad", CuTrapEta, fiberClad, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//
//			 		G4LogicalVolume *FiberClad_Logic_ith
//			 		= new G4LogicalVolume(FiberClad_ith, clad_C_Material, "FiberClad_Logic");
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberClad_Logic_ith, "FiberClad_Phys", CuTrapEta_Logical, false, TrapEta_i_total, checkOverlaps);
//
//			 		// -- Cores -- //
//			 		G4VSolid* FiberCore_ith;
//			 		G4LogicalVolume *FiberCore_Logic_ith;
//
//			 		// -- s, c, s, c, ... -- //
//			 		// -- c, s, c, s, ... -- //
//			 		bool isFiberC = this->IsFiberC(i_EtaDir, i_PhiDir);
//			 		if( isFiberC )
//			 		{
//			 			++TrapEta_fiberC_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuTrapEta, fiberCoreC, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_C_Material, "FiberCore_Logic");
//			 		}
//			 		else
//			 		{
//			 			++TrapEta_fiberS_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuTrapEta, fiberCoreS, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_S_Material, "FiberCore_Logic");
//			 		}
//
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberCore_Logic_ith, "FiberCore_Phys", FiberClad_Logic_ith, false, TrapEta_i_total, checkOverlaps);
//
//			 		G4VisAttributes* visAttr = new G4VisAttributes();
//			 		if( isFiberC )
//			 			visAttr->SetColour( G4Colour(0.0,0.0,1.0) ); // -- blue -- //
//			 		else
//			 			visAttr->SetColour( G4Colour(1.0,1.0,0.0) );  // -- yellow -- //
//			 		visAttr->SetForceSolid(true);
//			 		visAttr->SetVisibility(true);
//			 		FiberCore_Logic_ith->SetVisAttributes(visAttr);
//			 	}
//			 }
//			 cout<<i_barrel<<" CuTrapEta(rectangular part)"<<endl;
//			 cout<<"Cerenkov fiber #: "<<TrapEta_fiberC_total<<"		Scintillator fiber #: "<<TrapEta_fiberS_total<<"		total fiber #: "<<TrapEta_i_total<<endl;
//
//			//////////////////////////////////////////////////////////
//			//-- Cu trapezoidEta for triangular part fiber-- //
//			/////////////////////////////////////////////////////////
//
//
//			// -- Trdfibers for triangular part-- //
//
//			G4double CuTrapTrdLen_EtaDir = CuTrapLen_EtaDir*cos(CuTrdAngle_face);
//			const G4int TrapTrd_nFiber_EtaDir = floor( CuTrapTrdLen_EtaDir /dist_btwCore) - 1;
//			G4double TrapTrd_dist_edge_EtaDir = ( CuTrapTrdLen_EtaDir - (TrapTrd_nFiber_EtaDir-1)*dist_btwCore)/ 2.0;
//			const G4int TrapTrd_nFiber_PhiDirOuter = floor(CuTrdLen_PhiDirShort_Pre / dist_btwCore ) - 1;
//			G4double TrapTrd_dist_edge_PhiDirOuter = ( CuTrdLen_PhiDirShort_Pre - (TrapTrd_nFiber_PhiDirOuter-1)*dist_btwCore ) / 2.0;
//					
//			G4int TrapTrdEta_i_total = 0;
//			G4int TrapTrdEta_fiberC_total = 0;
//			G4int TrapTrdEta_fiberS_total = 0;
//			for(G4int i_EtaDir=0; i_EtaDir<TrapTrd_nFiber_EtaDir; i_EtaDir++)
//			{
//
//				G4double CuTrapTrdLen_PhiDir_Fiber
//				= CuTrdLen_PhiDirLong + 2*(TrapTrd_dist_edge_EtaDir + dist_btwCore*i_EtaDir)*tan(CuTrdAngle_face);
//				const G4int TrapTrd_nFiber_PhiDir = floor( CuTrapTrdLen_PhiDir_Fiber / dist_btwCore ) - 1;	
//				
//			 	for(G4int i_PhiDir=0; i_PhiDir<TrapTrd_nFiber_PhiDir; i_PhiDir++)
//			 	{
//			 		TrapTrdEta_i_total++;
//			 		G4double x_EtaDir = CuLen_EtaDir/2.0 + (TrapTrd_dist_edge_EtaDir + dist_btwCore*i_EtaDir)/cos(CuTrdAngle_face) - (TrapTrd_dist_edge_PhiDirOuter + dist_btwCore*i_PhiDir)*sin(CuTrdAngle_face);
//			 		G4double x_PhiDir = CuLen_PhiDir/2.0 + (TrapTrd_dist_edge_PhiDirOuter + dist_btwCore*i_PhiDir)*cos(CuTrdAngle_face);
//			 		// -- cladding: same shape for both C and S fiber -- //
//			 		G4VSolid* FiberClad_ith
//			 		= new G4IntersectionSolid("fiberClad", CuTrapTrdEta, fiberClad, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//
//			 		G4LogicalVolume *FiberClad_Logic_ith
//			 		= new G4LogicalVolume(FiberClad_ith, clad_C_Material, "FiberClad_Logic");
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberClad_Logic_ith, "FiberClad_Phys", CuTrapTrdEta_Logical, false, TrapTrdEta_i_total, checkOverlaps);
//
//			 		// -- Cores -- //
//			 		G4VSolid* FiberCore_ith;
//			 		G4LogicalVolume *FiberCore_Logic_ith;
//
//			 		// -- s, c, s, c, ... -- //
//			 		// -- c, s, c, s, ... -- //
//			 		bool isFiberC = this->IsFiberC(i_EtaDir, i_PhiDir);
//			 		if( isFiberC )
//			 		{
//			 			++TrapTrdEta_fiberC_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuTrapTrdEta, fiberCoreC, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_C_Material, "FiberCore_Logic");
//			 		}
//			 		else
//			 		{
//			 			++TrapTrdEta_fiberS_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuTrapTrdEta, fiberCoreS, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_S_Material, "FiberCore_Logic");
//			 		}
//
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberCore_Logic_ith, "FiberCore_Phys", FiberClad_Logic_ith, false, TrapTrdEta_i_total, checkOverlaps);
//
//			 		G4VisAttributes* visAttr = new G4VisAttributes();
//			 		if( isFiberC )
//			 			visAttr->SetColour( G4Colour(0.0,0.0,1.0) ); // -- blue -- //
//			 		else
//			 			visAttr->SetColour( G4Colour(1.0,1.0,0.0) );  // -- yellow -- //
//			 		visAttr->SetForceSolid(true);
//			 		visAttr->SetVisibility(true);
//			 		FiberCore_Logic_ith->SetVisAttributes(visAttr);
//			 	}
//			 }
//			cout<<i_barrel<<" CuTrapEta(for Phi trapezoid)"<<endl;
//			cout<<"Cerenkov fiber #: "<<TrapTrdEta_fiberC_total<<"		Scintillator fiber #: "<<TrapTrdEta_fiberS_total<<"		total fiber #: "<<TrapTrdEta_i_total<<endl;
//			fiberNumber += i_total + tri_i_total + TrapEta_i_total + TrapTrdEta_i_total;
//			cout<<fiberNumber<<endl;
//			cout<<endl;



			// G4ThreeVector origin(x,y,z);
			// G4RotationMatrix* RotMatrix = new G4RotationMatrix();

			// RotMatrix->rotateZ(90*deg);
			// RotMatrix->rotateZ(-i*theta_unit_ZRot);
			// RotMatrix->rotateX(90*deg);
			// RotMatrix->rotateX(-theta_unit*(copyNo+0.5));

			// -- place it -- //
			// new G4PVPlacement( RotMatrix, G4ThreeVector(), CuLogical, "CuPhysical", worldLogical, false, 0, checkOverlaps );
		}
	}

	///////////////////////////////
	// -- EndCap Geometry -- //
	//////////////////////////////

	G4double EndCap_CuLen_EtaDir =  CuLen_EtaDir;
	G4double EndCap_radius = sqrt(pow(2.5*m,2) + pow(radius,2));
	G4double end_dEta = 2*atan(EndCap_CuLen_EtaDir/(2.0*EndCap_radius));
	cout<<"cosine value : "<<cos(90*deg - Eta_Max - 40.5*end_dEta)<<endl;
	for(G4int i_cu=-41; i_cu<41; i_cu++)
	{
		G4double dEta;
		G4double dEta_Pre;
		if (i_cu==0) continue;
		else if (i_cu>0) {
			dEta = Eta_Max + (i_cu - 0.5)*end_dEta;
			dEta_Pre = Eta_Max + (i_cu - 0.5)*end_dEta;
		}
		else {
			dEta = -(Eta_Max - (i_cu + 0.5)*end_dEta);
			dEta_Pre = -(Eta_Max - (i_cu + 0.5)*end_dEta);
		}
		/////////////////////////////////////////
		// -- Logical Volume of Cu box -- //
		////////////////////////////////////////
		G4double Phi_radius_Top = EndCap_radius*cos(dEta) + EndCap_CuLen_EtaDir*abs(sin(dEta))/2.0;
		G4double Phi_radius_Bottom = EndCap_radius*cos(dEta) - EndCap_CuLen_EtaDir*abs(sin(dEta))/2.0;
		G4double Phi_radius_Pre = EndCap_radius*cos(dEta_Pre) - EndCap_CuLen_EtaDir*abs(sin(dEta_Pre))/2.0;
		G4double EndCap_dPhi = (2*pi) / nTower_PhiDir;
		G4double EndCap_half_dPhi = 0.5*EndCap_dPhi;
		G4double EndCap_CuLen_PhiDir_Top = 2*Phi_radius_Top*std::tan(EndCap_half_dPhi);
		G4double EndCap_CuLen_PhiDir_Bottom = 2*Phi_radius_Bottom*std::tan(EndCap_half_dPhi);
		G4double EndCap_CuLen_PhiDir_Pre = 2*Phi_radius_Pre*std::tan(EndCap_half_dPhi);

		vector<G4TwoVector> Box_Vertices;
		Box_Vertices.push_back(G4TwoVector(-EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Bottom/2.0));
		Box_Vertices.push_back(G4TwoVector(-EndCap_CuLen_EtaDir/2.0,-EndCap_CuLen_PhiDir_Bottom/2.0));
		Box_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,-EndCap_CuLen_PhiDir_Top/2.0));
		Box_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		Box_Vertices.push_back(G4TwoVector(-EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Bottom/2.0));
		Box_Vertices.push_back(G4TwoVector(-EndCap_CuLen_EtaDir/2.0,-EndCap_CuLen_PhiDir_Bottom/2.0));
		Box_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,-EndCap_CuLen_PhiDir_Top/2.0));
		Box_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		
		G4GenericTrap *EndCap_CuBox
		= new G4GenericTrap("EndCap_CuBox", CuLen_H/2.0,Box_Vertices);
		G4LogicalVolume *EndCap_CuLogical
		= new G4LogicalVolume(EndCap_CuBox, cu, "EndCap_CuLogical");
		
		//////////////////////////////////////////////////////////////////////////////////
		// -- Logical volume of Cu trapezoid for Eta direction -- //
		////////////////////////////////////////////////////////////////////////////////
		
		G4double EndCap_CuTrdLen_EtaDir = CuLen_H * tan(end_dEta);
		
		vector<G4TwoVector> Trd_Vertices;
		Trd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		Trd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,-EndCap_CuLen_PhiDir_Top/2.0));
		Trd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,-EndCap_CuLen_PhiDir_Top/2.0));
		Trd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		Trd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		Trd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,-EndCap_CuLen_PhiDir_Top/2.0));
		Trd_Vertices.push_back(G4TwoVector(EndCap_CuTrdLen_EtaDir + EndCap_CuLen_EtaDir/2.0,-EndCap_CuLen_PhiDir_Pre/2.0));
		Trd_Vertices.push_back(G4TwoVector(EndCap_CuTrdLen_EtaDir + EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Pre/2.0));
		G4GenericTrap *EndCap_CuTrd
		= new G4GenericTrap("EndCap_CuTrd", CuLen_H/2.0,Trd_Vertices);
		G4LogicalVolume *EndCap_CuTrdLogical
		= new G4LogicalVolume(EndCap_CuTrd, cu, "EndCap_CuTrdLogical");
		
		/////////////////////////////////////////////////////////////////////
		// -- Logical volume of Cu trapezoid for Phi direction -- //
		////////////////////////////////////////////////////////////////////
		
		G4double Angle_btwTower_cos = pow(cos(dEta),2)*cos(EndCap_dPhi) + pow(sin(dEta),2);
		G4double Angle_btwTower_tan = sqrt(1-pow(Angle_btwTower_cos,2))/Angle_btwTower_cos;
		G4double Box_diffOfPhiDir = EndCap_CuLen_PhiDir_Top - EndCap_CuLen_PhiDir_Bottom;
		G4double CuBoxAngle_face = atan(Box_diffOfPhiDir/(2.0*EndCap_CuLen_EtaDir));
		G4double CuTrapAngle_face = acos((pow(sin(dEta),2)*cos(EndCap_dPhi) + pow(cos(dEta),2))*pow(cos(CuBoxAngle_face),2) + pow(sin(CuBoxAngle_face),2)*cos(EndCap_dPhi))/2.0;
		G4double CuTrapLen_PhiDir = CuLen_H*Angle_btwTower_tan;
		
		vector<G4TwoVector> Trap_Vertices;
		Trap_Vertices.push_back(G4TwoVector(-EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Bottom/2.0));
		Trap_Vertices.push_back(G4TwoVector(-EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Bottom/2.0));
		Trap_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		Trap_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		Trap_Vertices.push_back(G4TwoVector(-EndCap_CuLen_EtaDir/2.0 - CuTrapLen_PhiDir*sin(CuTrapAngle_face),CuTrapLen_PhiDir*cos(CuTrapAngle_face) + EndCap_CuLen_PhiDir_Bottom/2.0));
		Trap_Vertices.push_back(G4TwoVector(-EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Bottom/2.0));
		Trap_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		Trap_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0 - CuTrapLen_PhiDir*sin(CuTrapAngle_face),CuTrapLen_PhiDir*cos(CuTrapAngle_face) + EndCap_CuLen_PhiDir_Top/2.0));
		G4GenericTrap *EndCap_CuTrap
		= new G4GenericTrap("EndCap_CuTrap", CuLen_H/2.0,Trap_Vertices);
		G4LogicalVolume *EndCap_CuTrapLogical
		= new G4LogicalVolume(EndCap_CuTrap, cu, "EndCap_CuTrapLogical");

		/////////////////////////////////////////////////////////////////////
		// -- Logical volume of Cu trapezoid for diagonal direction -- //
		////////////////////////////////////////////////////////////////////
		
		vector<G4TwoVector> TrapTrd_Vertices;
		TrapTrd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		TrapTrd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		TrapTrd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		TrapTrd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		TrapTrd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0 - CuTrapLen_PhiDir*sin(CuTrapAngle_face),CuTrapLen_PhiDir*cos(CuTrapAngle_face) + EndCap_CuLen_PhiDir_Top/2.0));
		TrapTrd_Vertices.push_back(G4TwoVector(EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Top/2.0));
		TrapTrd_Vertices.push_back(G4TwoVector(EndCap_CuTrdLen_EtaDir + EndCap_CuLen_EtaDir/2.0,EndCap_CuLen_PhiDir_Pre/2.0));
		TrapTrd_Vertices.push_back(G4TwoVector(EndCap_CuTrdLen_EtaDir*cos(2*CuTrapAngle_face) + EndCap_CuLen_EtaDir/2.0 - CuTrapLen_PhiDir*sin(CuTrapAngle_face),EndCap_CuTrdLen_EtaDir*sin(2*CuTrapAngle_face) + CuTrapLen_PhiDir*cos(CuTrapAngle_face) + EndCap_CuLen_PhiDir_Pre/2.0));
		
		G4GenericTrap *EndCap_CuTrapTrd
		= new G4GenericTrap("EndCap_CuTrapTrd", CuLen_H/2.0,TrapTrd_Vertices);
		G4LogicalVolume *EndCap_CuTrapTrdLogical
		= new G4LogicalVolume(EndCap_CuTrapTrd, cu, "EndCap_CuTrapTrdLogical");

		for (G4int i_Phi=0;i_Phi<2;i_Phi++){
			
			////////////////////////////////////
			// -- Placement of Cu box -- //
			///////////////////////////////////
			G4double CuLen_CenterH = EndCap_radius + CuLen_H/2.0;
			G4double dPhi = i_Phi * EndCap_dPhi;
			G4RotationMatrix rotM  = G4RotationMatrix();
			if (i_cu<0) rotM.rotateZ(180*deg);
			rotM.rotateY(90*deg -dEta);
			rotM.rotateZ(dPhi);
			G4double Phi_move = CuLen_CenterH*cos(dEta);
			G4double Eta_move = CuLen_CenterH*sin(dEta);
			G4ThreeVector Trans_Vector = G4ThreeVector(Phi_move*cos(dPhi), Phi_move*sin(dPhi) , Eta_move);
			G4Transform3D transform = G4Transform3D(rotM, Trans_Vector);
			new G4PVPlacement(transform, EndCap_CuLogical, "EndCap_CuPhysical", worldLogical, false, i_cu, checkOverlaps );
			
			////////////////////////////////////////////////////////////////////////////
			// -- Placement of Cu trapezoid for Eta direction (triangle) -- //
			///////////////////////////////////////////////////////////////////////////
			new G4PVPlacement(transform, EndCap_CuTrdLogical, "EndCap_CuTrdPhysical", worldLogical, false, i_cu, checkOverlaps );
			
			new G4PVPlacement(transform, EndCap_CuTrapLogical, "EndCap_CuTrapPhysical", worldLogical, false, i_cu, checkOverlaps );
			new G4PVPlacement(transform, EndCap_CuTrapTrdLogical, "EndCap_CuTrapTrdPhysical", worldLogical, false, i_cu, checkOverlaps );
			


//			////////////////////////////////////
//			// -- EndCap Cu box fiber -- //
//			///////////////////////////////////
//			 G4int i_total = 0;
//			 G4int fiberC_total = 0;
//			 G4int fiberS_total = 0;
//			 const G4int nFiber_EtaDir = floor( EndCap_CuLen_EtaDir / dist_btwCore ) - 1;
//			 G4double dist_edge_EtaDir = ( EndCap_CuLen_EtaDir - (nFiber_EtaDir-1)*dist_btwCore ) / 2.0;
//			 const G4int nFiber_PhiDirOuter = floor( EndCap_CuLen_PhiDir_Top / dist_btwCore ) - 1;
//			 G4double dist_edge_PhiDirOuter = ( EndCap_CuLen_PhiDir_Top - (nFiber_PhiDirOuter-1)*dist_btwCore ) / 2.0;
//
//			 for(G4int i_EtaDir=0; i_EtaDir<nFiber_EtaDir; i_EtaDir++)
//			 {
//				const G4int nFiber_PhiDir = floor( (EndCap_CuLen_PhiDir_Bottom + 2*(dist_edge_EtaDir + dist_btwCore*i_EtaDir)*tan(CuBoxAngle_face))/ dist_btwCore ) - 1;
//			 	G4double x_EtaDir = ((-1)*EndCap_CuLen_EtaDir / 2.0) + dist_edge_EtaDir + dist_btwCore*i_EtaDir;
//
//				for(G4int i_PhiDir=0; i_PhiDir<nFiber_PhiDir; i_PhiDir++)
//			 	{
//			 		i_total++;
//					G4double x_PhiDir = ((-1)*EndCap_CuLen_PhiDir_Bottom / 2.0) - (dist_edge_EtaDir + dist_btwCore*i_EtaDir)*tan(CuBoxAngle_face) + dist_edge_PhiDirOuter + dist_btwCore*i_PhiDir;
//
//			 		// -- cladding: same shape for both C and S fiber -- //
//			 		G4VSolid* FiberClad_ith
//			 		= new G4IntersectionSolid("fiberClad", EndCap_CuBox, fiberClad, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//
//			 		G4LogicalVolume *FiberClad_Logic_ith
//			 		= new G4LogicalVolume(FiberClad_ith, clad_C_Material, "FiberClad_Logic");
//
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberClad_Logic_ith, "FiberClad_Phys", EndCap_CuLogical, false, i_total, checkOverlaps);
//
//			 		// -- Cores -- //
//			 		G4VSolid* FiberCore_ith;
//			 		G4LogicalVolume *FiberCore_Logic_ith;
//			 		// -- s, c, s, c, ... -- //
//			 		// -- c, s, c, s, ... -- //
//			 		bool isFiberC = this->IsFiberC(i_EtaDir, i_PhiDir);
//			 		if( isFiberC )
//			 		{
//			 			++fiberC_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", EndCap_CuBox, fiberCoreC, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_C_Material, "FiberCore_Logic");
//			 		}
//			 		else
//			 		{
//			 			++fiberS_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", EndCap_CuBox, fiberCoreS, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_S_Material, "FiberCore_Logic");
//			 		}
//
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberCore_Logic_ith, "FiberCore_Phys", FiberClad_Logic_ith, false, i_total, checkOverlaps);
//
//			 		G4VisAttributes* visAttr = new G4VisAttributes();
//			 		if( isFiberC )
//			 			visAttr->SetColour( G4Colour(0.0,0.0,1.0) ); // -- blue -- //
//			 		else
//			 			visAttr->SetColour( G4Colour(1.0,1.0,0.0) );  // -- yellow -- //
//			 		visAttr->SetForceSolid(true);
//			 		visAttr->SetVisibility(true);
//			 		FiberCore_Logic_ith->SetVisAttributes(visAttr);
//			 	}
//			 }
//			 cout<<i_cu<<" CuBox"<<endl;
//			 cout<<"Cerenkov fiber #: "<<fiberC_total<<"		Scintillator fiber #: "<<fiberS_total<<"		total fiber #: "<<i_total<<endl;
//
//			//////////////////////////////////////////////////////////////
// 			// -- EndCap Cu Phi trapezoid fiber (called Trap)-- //
// 			/////////////////////////////////////////////////////////////
//
//			G4double CuTrdLen_EtaDir = CuLen_EtaDir*cos(CuTrapAngle_face);
// 			const G4int tri_nFiber_EtaDir = floor( CuTrdLen_EtaDir / dist_btwCore ) - 1;
// 			G4double tri_dist_edge_EtaDir = ( CuTrdLen_EtaDir - (tri_nFiber_EtaDir-1)*dist_btwCore ) / 2.0;
// 			const G4int tri_nFiber_PhiDir = floor( CuTrapLen_PhiDir / dist_btwCore ) - 1;
// 			G4double tri_dist_edge_PhiDir = ( CuTrapLen_PhiDir - (tri_nFiber_PhiDir-1)*dist_btwCore ) / 2.0;
//
// 			G4int tri_i_total = 0;
// 			G4int tri_fiberC_total = 0;
// 			G4int tri_fiberS_total = 0;
// 			for(G4int i_EtaDir=0; i_EtaDir< tri_nFiber_EtaDir; i_EtaDir++)
// 			{
// 				 for(G4int i_PhiDir=0; i_PhiDir<tri_nFiber_PhiDir; i_PhiDir++)
// 				 {
// 				 	tri_i_total++;
// 				 	G4double x_EtaDir
// 				 	= ((-1)*CuLen_EtaDir / 2.0) + (tri_dist_edge_EtaDir + dist_btwCore*i_EtaDir)*cos(CuTrapAngle_face) + (tri_dist_edge_PhiDir + dist_btwCore*i_PhiDir)*sin(CuTrapAngle_face);
// 				 	G4double x_PhiDir
// 					= EndCap_CuLen_PhiDir_Bottom/2.0 + (tri_dist_edge_EtaDir + dist_btwCore*i_EtaDir)*sin(CuTrapAngle_face) + (tri_dist_edge_PhiDir + dist_btwCore*i_PhiDir)*cos(CuTrapAngle_face);
//
// 				 	G4ThreeVector trans_vector=G4ThreeVector(x_EtaDir, x_PhiDir, 0);
//
// 				 	// -- cladding: same shape for both C and S fiber -- //
// 			 		G4VSolid *FiberClad_ith
// 					= new G4IntersectionSolid("fiberClad", EndCap_CuTrap, fiberClad, 0, trans_vector);
//
// 			 		G4LogicalVolume *FiberClad_Logic_ith
// 			 		= new G4LogicalVolume(FiberClad_ith, clad_C_Material, "FiberClad_Logic");
// 			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberClad_Logic_ith, "FiberClad_Phys", EndCap_CuTrapLogical, false, tri_i_total, checkOverlaps);
//
//			 		// -- Cores -- //
// 			 		G4VSolid* FiberCore_ith;
// 			 		G4LogicalVolume *FiberCore_Logic_ith;
// 				 	// -- s, c, s, c, ... -- //
// 			 		// -- c, s, c, s, ... -- //
// 			 		bool isFiberC = this->IsFiberC(i_EtaDir, i_PhiDir);
// 			 		if( isFiberC )
// 			 		{
// 			 			++tri_fiberC_total;
// 			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", EndCap_CuTrap, fiberCoreC, 0, trans_vector);
// 			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_C_Material, "FiberCore_Logic");
// 			 		}
// 			 		else
// 			 		{
// 			 			++tri_fiberS_total;
// 			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", EndCap_CuTrap, fiberCoreS, 0, trans_vector);
// 			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_S_Material, "FiberCore_Logic");
// 			 		}
// 			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberCore_Logic_ith, "FiberCore_Phys", FiberClad_Logic_ith, false, tri_i_total, checkOverlaps);
// 				 		G4VisAttributes* visAttr1 = new G4VisAttributes();
// 			 		if( isFiberC ){
// 			 			visAttr1->SetColour( G4Colour(0.0,0.0,1.0) ); // -- blue -- //
// 			 		}
// 			 		else{
// 			 			visAttr1->SetColour( G4Colour(1.0,1.0,0.0) );  // -- yellow -- //
// 			 		}
// 			 		visAttr1->SetForceSolid(true);
// 			 		visAttr1->SetVisibility(true);
// 			 		FiberCore_Logic_ith->SetVisAttributes(visAttr1);
// 			 	}
// 			 }
// 		cout<<i_cu<<" CuTrapPhi"<<endl;
// 		cout<<"Cerenkov fiber #: "<<tri_fiberC_total<<"		Scintillator fiber #: "<<tri_fiberS_total<<"		total fiber #: "<<tri_i_total<<endl;
//
// 			////////////////////////////////////////////////////////////
//			//-- EndCap Cu Eta trapezoid fiber (called Trd)-- //
//			//////////////////////////////////////////////////////////
//
//			const G4int Trd_nFiber_EtaDir = floor( EndCap_CuTrdLen_EtaDir / dist_btwCore) - 1;
//			G4double Trd_dist_edge_EtaDir = ( EndCap_CuTrdLen_EtaDir - (Trd_nFiber_EtaDir-1)*dist_btwCore)/ 2.0;
//			const G4int Trd_nFiber_PhiDirOuter = floor( EndCap_CuLen_PhiDir_Pre / dist_btwCore ) - 1;
//			G4double Trd_dist_edge_PhiDirOuter = ( EndCap_CuLen_PhiDir_Pre - (Trd_nFiber_PhiDirOuter-1)*dist_btwCore ) / 2.0;
//			G4int TrapEta_i_total = 0;
//			G4int TrapEta_fiberC_total = 0;
//			G4int TrapEta_fiberS_total = 0;
//
//			for(G4int i_EtaDir=0; i_EtaDir<Trd_nFiber_EtaDir; i_EtaDir++)
//			{
//				const G4int Trd_nFiber_PhiDir = floor(( EndCap_CuLen_PhiDir_Top + 2*(Trd_dist_edge_EtaDir + dist_btwCore*i_EtaDir)*tan(CuBoxAngle_face)) / dist_btwCore ) - 1;
//				G4double x_EtaDir = EndCap_CuLen_EtaDir / 2.0 + Trd_dist_edge_EtaDir + dist_btwCore*i_EtaDir;
//			 	for(G4int i_PhiDir=0; i_PhiDir<Trd_nFiber_PhiDirOuter; i_PhiDir++)
//			 	{
//			 		TrapEta_i_total++;
//			 		G4double x_PhiDir = ((-1)*EndCap_CuLen_PhiDir_Top / 2.0)  - (Trd_dist_edge_EtaDir + dist_btwCore*i_EtaDir)*tan(CuBoxAngle_face) + Trd_dist_edge_PhiDirOuter + dist_btwCore*i_PhiDir;
//
//			 		// -- cladding: same shape for both C and S fiber -- //
//			 		G4VSolid* FiberClad_ith
//			 		= new G4IntersectionSolid("fiberClad", EndCap_CuTrd, fiberClad, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//
//			 		G4LogicalVolume *FiberClad_Logic_ith
//			 		= new G4LogicalVolume(FiberClad_ith, clad_C_Material, "FiberClad_Logic");
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberClad_Logic_ith, "FiberClad_Phys", EndCap_CuTrdLogical, false, TrapEta_i_total, checkOverlaps);
//
//			 		// -- Cores -- //
//			 		G4VSolid* FiberCore_ith;
//			 		G4LogicalVolume *FiberCore_Logic_ith;
//
//			 		// -- s, c, s, c, ... -- //
//			 		// -- c, s, c, s, ... -- //
//			 		bool isFiberC = this->IsFiberC(i_EtaDir, i_PhiDir);
//			 		if( isFiberC )
//			 		{
//			 			++TrapEta_fiberC_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", EndCap_CuTrd, fiberCoreC, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_C_Material, "FiberCore_Logic");
//			 		}
//			 		else
//			 		{
//			 			++TrapEta_fiberS_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", EndCap_CuTrd, fiberCoreS, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_S_Material, "FiberCore_Logic");
//			 		}
//
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberCore_Logic_ith, "FiberCore_Phys", FiberClad_Logic_ith, false, TrapEta_i_total, checkOverlaps);
//
//			 		G4VisAttributes* visAttr = new G4VisAttributes();
//			 		if( isFiberC )
//			 			visAttr->SetColour( G4Colour(0.0,0.0,1.0) ); // -- blue -- //
//			 		else
//			 			visAttr->SetColour( G4Colour(1.0,1.0,0.0) );  // -- yellow -- //
//			 		visAttr->SetForceSolid(true);
//			 		visAttr->SetVisibility(true);
//			 		FiberCore_Logic_ith->SetVisAttributes(visAttr);
//			 	}
//			 }
//			 cout<<i_cu<<" CuTrapEta(rectangular part)"<<endl;
//			 cout<<"Cerenkov fiber #: "<<TrapEta_fiberC_total<<"		Scintillator fiber #: "<<TrapEta_fiberS_total<<"		total fiber #: "<<TrapEta_i_total<<endl;
//			//////////////////////////////////////////////////////////////
//			//-- Cu diagonal trapezoid fiber (called TrapTrd)-- //
//			/////////////////////////////////////////////////////////////
//
//
//			// -- Trdfibers for triangular part-- //
//
//			G4double CuTrapTrdLen_EtaDir = EndCap_CuTrdLen_EtaDir*cos(CuTrapAngle_face);
//			const G4int TrapTrd_nFiber_EtaDir = floor( CuTrapTrdLen_EtaDir /dist_btwCore) - 1;
//			G4double TrapTrd_dist_edge_EtaDir = ( CuTrapTrdLen_EtaDir - (TrapTrd_nFiber_EtaDir-1)*dist_btwCore)/ 2.0;
//			const G4int TrapTrd_nFiber_PhiDir = floor(CuTrapLen_PhiDir / dist_btwCore ) - 1;
//			G4double TrapTrd_dist_edge_PhiDir = ( CuTrapLen_PhiDir - (TrapTrd_nFiber_PhiDir-1)*dist_btwCore ) / 2.0;
//
//			G4int TrapTrdEta_i_total = 0;
//			G4int TrapTrdEta_fiberC_total = 0;
//			G4int TrapTrdEta_fiberS_total = 0;
//			for(G4int i_EtaDir=0; i_EtaDir<TrapTrd_nFiber_EtaDir; i_EtaDir++)
//			{
//			 	for(G4int i_PhiDir=0; i_PhiDir<TrapTrd_nFiber_PhiDir; i_PhiDir++)
//			 	{
//			 		TrapTrdEta_i_total++;
//			 		G4double x_EtaDir = EndCap_CuLen_EtaDir/2.0 + (TrapTrd_dist_edge_EtaDir + dist_btwCore*i_EtaDir)*cos(CuTrapAngle_face) - (TrapTrd_dist_edge_PhiDir + dist_btwCore*i_PhiDir)*sin(CuTrapAngle_face);
//			 		G4double x_PhiDir = EndCap_CuLen_PhiDir_Top/2.0 + (TrapTrd_dist_edge_EtaDir + dist_btwCore*i_EtaDir)*sin(CuTrapAngle_face) + (TrapTrd_dist_edge_PhiDir + dist_btwCore*i_PhiDir)*cos(CuTrapAngle_face);
//			 		// -- cladding: same shape for both C and S fiber -- //
//			 		G4VSolid* FiberClad_ith
//			 		= new G4IntersectionSolid("fiberClad", EndCap_CuTrapTrd, fiberClad, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//
//			 		G4LogicalVolume *FiberClad_Logic_ith
//			 		= new G4LogicalVolume(FiberClad_ith, clad_C_Material, "FiberClad_Logic");
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberClad_Logic_ith, "FiberClad_Phys", EndCap_CuTrapTrdLogical, false, TrapTrdEta_i_total, checkOverlaps);
//
//			 		// -- Cores -- //
//			 		G4VSolid* FiberCore_ith;
//			 		G4LogicalVolume *FiberCore_Logic_ith;
//
//			 		// -- s, c, s, c, ... -- //
//			 		// -- c, s, c, s, ... -- //
//			 		bool isFiberC = this->IsFiberC(i_EtaDir, i_PhiDir);
//			 		if( isFiberC )
//			 		{
//			 			++TrapTrdEta_fiberC_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", EndCap_CuTrapTrd, fiberCoreC, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_C_Material, "FiberCore_Logic");
//			 		}
//			 		else
//			 		{
//			 			++TrapTrdEta_fiberS_total;
//			 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", EndCap_CuTrapTrd, fiberCoreS, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
//			 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_S_Material, "FiberCore_Logic");
//			 		}
//
//			 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberCore_Logic_ith, "FiberCore_Phys", FiberClad_Logic_ith, false, TrapTrdEta_i_total, checkOverlaps);
//
//			 		G4VisAttributes* visAttr = new G4VisAttributes();
//			 		if( isFiberC )
//			 			visAttr->SetColour( G4Colour(0.0,0.0,1.0) ); // -- blue -- //
//			 		else
//			 			visAttr->SetColour( G4Colour(1.0,1.0,0.0) );  // -- yellow -- //
//			 		visAttr->SetForceSolid(true);
//			 		visAttr->SetVisibility(true);
//			 		FiberCore_Logic_ith->SetVisAttributes(visAttr);
//			 	}
//			 }
//			cout<<i_cu<<" CuTrapEta(for Phi trapezoid)"<<endl;
//			cout<<"Cerenkov fiber #: "<<TrapTrdEta_fiberC_total<<"		Scintillator fiber #: "<<TrapTrdEta_fiberS_total<<"		total fiber #: "<<TrapTrdEta_i_total<<endl;
//			fiberNumber += i_total + tri_i_total + TrapEta_i_total + TrapTrdEta_i_total;
//			cout<<fiberNumber<<endl;
//			cout<<endl;
			}

		G4VisAttributes* visAttr1;
		visAttr1 = new G4VisAttributes( G4Colour(0.0,1.0,1.0) ); // -- cyan -- //
		visAttr1->SetVisibility(true);
		EndCap_CuLogical->SetVisAttributes(visAttr1);
		fVisAttributes.push_back(visAttr1);
	}
	// -- visualization -- //
	G4VisAttributes* visAttr1;
	visAttr1 = new G4VisAttributes( G4Colour(0.0,1.0,1.0) ); // -- cyan -- //
	visAttr1->SetVisibility(true);
	CuLogical->SetVisAttributes(visAttr1);
	fVisAttributes.push_back(visAttr1);
	G4VisAttributes* visAttr2 = new G4VisAttributes();
	visAttr2->SetColour( G4Colour(1.0,0.0,0.0) ); // -- red -- //
	visAttr2->SetForceSolid(true);
	visAttr2->SetVisibility(true);
	this->PMTPCBox_Logic->SetVisAttributes(visAttr2);
	fVisAttributes.push_back(visAttr2);
	G4VisAttributes* visAttr3 = new G4VisAttributes();
	visAttr3->SetColour( G4Colour(0.0,1.0,0.0) ); // -- green -- //
	// visAttr3->SetForceSolid(true);
	visAttr3->SetVisibility(true);
	PMTGlassBox_Logic->SetVisAttributes(visAttr3);
	fVisAttributes.push_back(visAttr3);
	
//	G4RotationMatrix rotM  = G4RotationMatrix();
//	rotM.rotateZ(dPhi);
//	G4Box* testBox1
//	= new G4Box("CuBox", CuLen_EtaDir/2.0, CuLen_PhiDir/2.0, CuLen_H/2.0);
//	G4Box* testBox2_new
//	= new G4Box("CuBox", CuLen_EtaDir/2.0, CuLen_PhiDir/2.0, CuLen_H/2.0);
//	G4Transform3D transform = G4Transform3D(rotM, G4ThreeVector(0,0,0));
//	G4Transform3D transform1 = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,radius + CuLen_H/2.0));
//	G4VSolid* testBox2
//	= new G4intersectionSolid("CuBox", worldSolid, testBox2_new, transform);
//	G4LogicalVolume *testBoxLogical1
//	= new G4LogicalVolume(testBox1, cu, "CuLogical");
//	G4LogicalVolume *testBoxLogical2
//	= new G4LogicalVolume(testBox2, cu, "CuLogical");
//	new G4PVPlacement(transform1, testBoxLogical1, "EndCap_CuPhysical", worldLogical, false, 0, checkOverlaps );
//	new G4PVPlacement(transform1, testBoxLogical2, "EndCap_CuPhysical", worldLogical, false, 0, checkOverlaps );
	
	
	return worldPhysical;
}

void WGR16DetectorConstruction::ConstructSDandField()
{
	G4SDManager* SDManager = G4SDManager::GetSDMpointer();

	G4String PMTName = "WGR16/PMTSD";
	G4String HitCollectionName = "PMTColl";
	G4VSensitiveDetector* SD_PMTPC = new WGR16PMTSD(PMTName, HitCollectionName);

	SDManager->AddNewDetector(SD_PMTPC);
	this->PMTPCBox_Logic->SetSensitiveDetector(SD_PMTPC);
}

void WGR16DetectorConstruction::ConstructMaterials()
{
	G4NistManager* nistManager = G4NistManager::Instance();
	
	// -- Vacuum: "Galactic" -- //
	nistManager->FindOrBuildMaterial("G4_Galactic");
	nistManager->FindOrBuildMaterial("G4_Cu");

	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

void WGR16DetectorConstruction::DefineCommands()
{
	// Define /WGR16/detector command directory using generic messenger class
	// fMessenger = new G4GenericMessenger(this, 
	//                                     "/WGR16/detector/", 
	//                                     "Detector control");

	// G4GenericMessenger::Command& B2ECollAngleCmd
	//   = fMessenger->DeclareMethodWithUnit("Block2Angle","deg",
	//       &WGR16DetectorConstruction::SetBlock2Angle,
	//       "Set rotation angle of exit collimator and corresponding detector");
	// B2ECollAngleCmd.SetParameterName("angle",true);
	// B2ECollAngleCmd.SetRange("angle>=-90. && angle<=90");
	// B2ECollAngleCmd.SetDefaultValue("60.");	
}

G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_PMMA()
{
	G4double PhotonEnergy[] = {
	2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

	// -- PMMA -- //
	G4double RefractiveIndex_PMMA[nEntries] =
	{
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
	};

	// G4double ABSLength_PMMA[nEntries] =
	// {
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m
	// };

	// G4double Reflectivity_PMMA[nEntries] =
	// {
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
	// };

	G4MaterialPropertiesTable* mpPMMA = new G4MaterialPropertiesTable();
	mpPMMA->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_PMMA,nEntries);
	//mpPMMA->AddProperty("ABSLENGTH",PhotonEnergy,ABSLength_PMMA,nEntries);
	//mpPMMA->AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity_PMMA,nEntries);

	return mpPMMA;
}


G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_FS()
{
	G4double PhotonEnergy[] = {
	2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

	//--- Fluorinated Polymer (FS) ---
	G4double RefractiveIndex_FluorinatedPolymer[nEntries] =
	{
	    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
	    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
	    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
	    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
	    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42
	};

	// G4double ABSLength_FS[nEntries] =
	// {
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m
	// };

	// G4double Reflectivity_FS[nEntries] =
	// {
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
	// };

	G4MaterialPropertiesTable* mpFS = new G4MaterialPropertiesTable();
	mpFS->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_FluorinatedPolymer,nEntries);

	//mpFS->AddProperty("ABSLENGTH",PhotonEnergy,ABSLength_FS,nEntries);
	//mpFS->AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity_FS,nEntries);

	return mpFS;
}

G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_Glass()
{
	G4double PhotonEnergy[] = {
	2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

	//
	//Glass
	//
	G4double RefractiveIndex_Glass[nEntries] =
	{   1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
	};
	
	// -- For test -- //
	// G4double RefractiveIndex_Glass[nEntries] =
	// {   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	//     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	//     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	//     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	//     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
	// };

	//Whenever I add some more detail optical properties, the GEANT4 will take lots of time to calculate
	//additional optical properties. I use only refractive indices to save time.
	// G4double Glass_AbsLength[nEntries] =
	// {
	//     420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,
	//     420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,
	//     420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,
	//     420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,
	//     420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,
	// };

	G4MaterialPropertiesTable* mpGlass = new G4MaterialPropertiesTable();
	mpGlass->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_Glass,nEntries);
	//mpGlass->AddProperty("ABSLENGTH",PhotonEnergy,Glass_AbsLength,nEntries);

	return mpGlass;
}

G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_Air()
{
	G4double PhotonEnergy[] = {
	2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

	//
	//Air
	//
	G4double RefractiveIndex_Air[nEntries] =
	{
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
	};

	G4MaterialPropertiesTable* mpAir = new G4MaterialPropertiesTable();
	mpAir->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Air, nEntries);

	return mpAir;
}

G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_PMTPC()
{
	G4double PhotonEnergy[] = {
	2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

	//
	//Air
	//
	G4double RefractiveIndex_Air[nEntries] =
	{
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
	};

	G4double p_mppc[2] = {2.00*eV, 3.47*eV};
	G4double refl_mppc[2] = {0.0, 0.0};
	G4double effi_mppc[2] = {0.11, 0.11}; // mimic Quantum Efficiency
	// G4double photocath_ReR[2] = {1.92, 1.92};
	// G4double photocath_ImR[2] = {1.69, 1.69};

	G4MaterialPropertiesTable* mpPMTPC = new G4MaterialPropertiesTable();
	mpPMTPC->AddProperty("REFLECTIVITY",p_mppc,refl_mppc,2);
	mpPMTPC->AddProperty("EFFICIENCY",p_mppc,effi_mppc,2);
	// mpPMTPC->AddProperty("REALINDEX",p_mppc,photocath_ReR,2);
	// mpPMTPC->AddProperty("IMAGINARYINDEX",p_mppc,photocath_ImR,2);
	mpPMTPC->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Air, nEntries);

	return mpPMTPC;
}

G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_PMTHouse()
{
	G4double ephoton[] = {2.00*eV, 3.47*eV};
	const G4int num = sizeof(ephoton) / sizeof(G4double);

	G4double reficy[] = {0.0,0.0};
	G4double efficy[] = {1.0,1.0};

	G4MaterialPropertiesTable* mpPMTHouse = new G4MaterialPropertiesTable();
	mpPMTHouse->AddProperty("REFLECTIVITY",ephoton,reficy,num);
	mpPMTHouse->AddProperty("EFFICIENCY",ephoton,efficy,num);

	return mpPMTHouse;
}


bool WGR16DetectorConstruction::IsFiberC(G4int i_EtaDir, G4int i_PhiDir)
{
	bool Flag = false;
	if( i_EtaDir % 2 == 0 ) // start with c fiber -- //
	{
		if( i_PhiDir % 2 == 0 ) // -- c fiber -- //
			Flag = true;
		else
			Flag = false;
	}
	else // -- start with s fiber -- //
	{
		if( i_PhiDir % 2 == 0 ) // -- c fiber -- //
			Flag = false;
		else
			Flag = true;

	}

	return Flag;
}
// Eta equation using proper geometry
G4double WGR16DetectorConstruction::Eta_Eq(G4double Pre_Eta, G4double Curr_Eta, G4double radius, G4double CuLen_EtaDir){
                return radius*(tan(Curr_Eta)-tan(Pre_Eta))+(CuLen_EtaDir/2)*(1/cos(Curr_Eta)-1/cos(Pre_Eta)-2*cos(Curr_Eta)-2*sin(Curr_Eta)*tan(Pre_Eta));
        }

// solving increasing equation by binary search
G4double WGR16DetectorConstruction::Solve_Eq(G4double Pre_Eta,G4double Low, G4double Max, G4double radius, G4double CuLen_EtaDir){
		if (Eta_Eq(Pre_Eta,Low,radius,CuLen_EtaDir)*Eta_Eq(Pre_Eta,Max,radius,CuLen_EtaDir)>=0){
			cout<<"no proper solution for Eta";
			return -1;
		}
       G4double eta;
       int cnt = 0;
       while (Max-Low>=0.0 && cnt<50){
    	   eta=(Low+Max)/2;
    	   G4double temp = Eta_Eq(Pre_Eta,eta,radius,CuLen_EtaDir);
    	   if (temp==0.0) return eta;
    	   else if (temp*Eta_Eq(Pre_Eta,Low,radius,CuLen_EtaDir)<0) Max = eta;
    	   else Low = eta;
    	   ++cnt;
       }
      return eta;
}



