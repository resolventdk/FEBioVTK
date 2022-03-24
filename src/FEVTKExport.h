#pragma once
#include <FECore/FECoreBase.h>
#include <FECore/FEModel.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>

class FEVTKExport
{
public:
	FEVTKExport(FEModel*);
	~FEVTKExport(void);

	bool Save();

private:
		
	bool AddDomains(vtkSmartPointer<vtkMultiBlockDataSet>);
	bool AddSurfaces(vtkSmartPointer<vtkMultiBlockDataSet>);

private:
	FEModel* m_fem;
	string m_prefix = "output";
	int m_ndump = 0;

	// check available names in FEBioMech/FEBioMechModule.cpp
	list<string> cell_data_fields {"stress", "Lagrange strain"};
	list<string> point_data_fields{};
	list<string> surface_data_fields{"contact traction", "contact gap"};

};