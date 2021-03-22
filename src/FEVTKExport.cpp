/*This file is part of the FEBio Studio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio-Studio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include "FEVTKExport.h"
#include <stdio.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEAnalysis.h>
#include <FEBioPlot/FEBioPlotFile.h>
#include "FECore/FECore.h"
#include <FECore/log.h>

enum VTK_CELLTYPE {
	VTK_VERTEX =                 1,
	VTK_POLY_VERTEX =            2,
	VTK_LINE =                   3,
	VTK_POLY_LINE =              4,
	VTK_TRIANGLE =               5,
	VTK_TRIANGLE_STRIP =         6,
	VTK_POLYGON =                7,
	VTK_PIXEL =                  8,
	VTK_QUAD =                   9,
	VTK_TETRA =                  10,
	VTK_VOXEL =                  11,
	VTK_HEXAHEDRON =             12,
	VTK_WEDGE =                  13,
	VTK_PYRAMID =                14,
	VTK_QUADRATIC_EDGE =         21,
	VTK_QUADRATIC_TRIANGLE =     22,
	VTK_QUADRATIC_QUAD =         23,
	VTK_QUADRATIC_TETRA =        24,
	VTK_QUADRATIC_HEXAHEDRON =   25,
	VTK_QUADRATIC_WEDGE =        26
};

void Space2_(char* szname)
{
	int n = (int)strlen(szname);
	for (int i = 0; i<n; ++i)
		if (szname[i] == ' ') szname[i] = '_';
}


FEVTKExport::FEVTKExport(FEModel* fem)
{
	m_fem = fem;

	// initialize json based meta file following Paraview format
	// https://gitlab.kitware.com/paraview/paraview/blob/v5.5.0/Documentation/release/ParaView-5.5.0.md#json-based-new-meta-file-format-for-series-added
	m_fp = fopen("test.vtk.series", "wt");
	fprintf(m_fp, "{\n");
	fprintf(m_fp, "  \"file-series-version\" : \"1.0\",\n");
	fprintf(m_fp, "  \"files\" : [\n");
	fclose(m_fp);
	m_fp = nullptr;

}

FEVTKExport::~FEVTKExport(void)
{
	// finalize json based meta file
	m_fp = fopen("test.vtk.series", "a");
	fprintf(m_fp, "  ]\n");
	fprintf(m_fp, "}\n");
	fclose(m_fp);
	m_fp = nullptr;
}

bool FEVTKExport::Save()
{

	// get current step and time
	int step = m_fem->GetCurrentStep()->m_ntimesteps;
	double time = m_fem->GetTime().currentTime;

	// set vtk file name
	char szfile[255];
	sprintf(szfile, "test%3.3d.vtk", step);

	// add this step to the meta file
	m_fp = fopen("test.vtk.series", "a");
	fprintf(m_fp, "     { \"name\" : \"%s\", \"time\" : %g },\n", szfile, time);
	fclose(m_fp);

	m_fp = fopen(szfile, "wt");
	if (m_fp == 0) return false;
        
	// --- H E A D E R ---
	WriteHeader();
	
	// --- N O D E S ---
	WritePoints();
        
    // --- E L E M E N T S ---
	WriteCells();
        
    // --- N O D E   D A T A ---
	printf("Write point data.\n");
	WritePointData();

    // --- E L E M E N T   C E L L   D A T A ---
	printf("Write cell data.\n");
	WriteCellData();
        
    fclose(m_fp);
	m_fp = nullptr;

	return true;
}

//-----------------------------------------------------------------------------
void FEVTKExport::WriteHeader()
{
	fprintf(m_fp, "%s\n"       ,"# vtk DataFile Version 3.0");
	fprintf(m_fp, "%s %g\n"    ,"vtk output at time", m_fem->GetCurrentTime());
	fprintf(m_fp, "%s\n"       ,"ASCII");
	fprintf(m_fp, "%s\n"       ,"DATASET UNSTRUCTURED_GRID");
}

//-----------------------------------------------------------------------------
void FEVTKExport::WritePoints()
{
	FEMesh& mesh = m_fem->GetMesh();
	int nodes = mesh.Nodes();
	fprintf(m_fp, "POINTS %d float\n", nodes);
	for (int j=0; j<nodes; j += 3)
	{
	    for (int k =0; k<3 && j+k<nodes;k++)
	    {
	        vec3d& r = mesh.Node(j+k).m_r0;
	        fprintf(m_fp, "%g %g %g ", r.x, r.y, r.z);
	    }
	    fprintf(m_fp, "\n");
	}
	fprintf(m_fp, "%s\n" ,"");
}

//-----------------------------------------------------------------------------
void FEVTKExport::WriteCells()
{
	FEMesh& mesh = m_fem->GetMesh();

	// get number of elements and node connections
	int NE = mesh.Elements();
	int nsize = 0;
	for (int j = 0; j < mesh.Elements(); ++j)
        nsize += mesh.Element(j)->Nodes() + 1;
	
	//// loop "surfaces" and add number of surface elements and node connections
	//for (int j = 0; j < mesh.Surfaces(); ++j) {
	//	FESurface& Sj = mesh.Surface(j);
	//	NE += Sj.Elements();
	//	for (int i = 0; i < Sj.Elements(); ++i)
	//		nsize += Sj.Element(i)->Nodes() + 1;
	//}

	// Write CELLS
    fprintf(m_fp, "CELLS %d %d\n", NE, nsize);
	// volume elements
	for (int j=0; j < mesh.Elements(); ++j) 
    {
		FEElement* el = mesh.Element(j);
        fprintf(m_fp, "%d ", el->Nodes());
        for (int k=0; k<el->Nodes(); ++k) fprintf(m_fp, "%d ", el->m_node[k]);
        fprintf(m_fp, "\n");
    }
	//// surface elements
	//for (int j = 0; j < mesh.Surfaces(); ++j) { // note surfaces may be dublicated
	//	FESurface& Sj = mesh.Surface(j);
	//	for (int i = 0; i < Sj.Elements(); ++i) {
	//		FEElement* el = Sj.Element(i);
	//		fprintf(m_fp, "%d ", el->Nodes());
	//		for (int k = 0; k < el->Nodes(); ++k) fprintf(m_fp, "%d ", el->m_node[k]);
	//		fprintf(m_fp, "\n");
	//	}
	//}
     
	// Write CELL_TYPES
    fprintf(m_fp, "\nCELL_TYPES %d\n", NE);
	int vtk_type;
	// volume elements
	for (int j = 0; j<mesh.Elements(); ++j)
    {
		FEElement* el = mesh.Element(j);
        switch (el->Type()) {
			//FE_HEX8G8,
			//FE_HEX8RI,
			//FE_HEX8G1,
			//FE_HEX20G8,
			//FE_HEX20G27,
			//FE_HEX27G27,
			//FE_TET4G1,
			//FE_TET4G4,
			//FE_TET5G4,
			//FE_TET10G1,
			//FE_TET10G4,
			//FE_TET10G8,
			//FE_TET10GL11,
			//FE_TET10G4RI1,
			//FE_TET10G8RI4,
			//FE_TET15G4,
			//FE_TET15G8,
			//FE_TET15G11,
			//FE_TET15G15,
			//FE_TET15G15RI4,
			//FE_TET20G15,
			//FE_PENTA6G6,
			//FE_PENTA15G8,
			//FE_PENTA15G21,
			//FE_PYRA5G8,
			//FE_PYRA13G8,
		case FE_TET4G1: vtk_type = VTK_TETRA; break;
		case FE_TET4G4: vtk_type = VTK_TETRA; break;
		default: printf("Unknown/unsupported volume element type!"); vtk_type = -1;  break;
        }            
        fprintf(m_fp, "%d\n", vtk_type);
    }
	//// surface elements
	//for (int j = 0; j < mesh.Surfaces(); ++j) {
	//	FESurface& Sj = mesh.Surface(j);
	//	for (int i = 0; i < Sj.Elements(); ++i) {
	//		FEElement* el = Sj.Element(i);
	//		switch (el->Type()) {
	//			//FE_QUAD4G4,
	//			//FE_QUAD4NI,
	//			//FE_TRI3G1,
	//			//FE_TRI3G3,
	//			//FE_TRI3G7,
	//			//FE_TRI3NI,
	//			//FE_TRI6G3,
	//			//FE_TRI6G4,
	//			//FE_TRI6G7,
	//			////	FE_TRI6MG7,
	//			//FE_TRI6GL7,
	//			//FE_TRI6NI,
	//			//FE_TRI7G3,
	//			//FE_TRI7G4,
	//			//FE_TRI7G7,
	//			//FE_TRI7GL7,
	//			//FE_TRI10G7,
	//			//FE_TRI10G12,
	//			//FE_QUAD8G9,
	//			//FE_QUAD8NI,
	//			//FE_QUAD9G9,
	//			//FE_QUAD9NI,
	//			case FE_TRI3G1 : vtk_type = VTK_TRIANGLE; break;
	//			case FE_TRI3G3 : vtk_type = VTK_TRIANGLE; break;
	//			case FE_TRI3G7 : vtk_type = VTK_TRIANGLE; break;
	//			case FE_TRI3NI : vtk_type = VTK_TRIANGLE; break;
	//		default: printf("Unknown/unsupported surface element type!"); vtk_type = -1;  break;
	//		}
	//		fprintf(m_fp, "%d\n", vtk_type);
	//	}
	//}
}

//-----------------------------------------------------------------------------
void FEVTKExport::WritePointData()
{
	FEMesh& mesh = m_fem->GetMesh();
	int N = mesh.Nodes();

		fprintf(m_fp, "\nPOINT_DATA %d\n", N);

		// write displacement - not a FEPlotData object
		fprintf(m_fp, "VECTORS displacement float\n");
		for (int i = 0; i < N; i++) {
			vec3d dr = mesh.Node(i).m_rt - mesh.Node(i).m_r0;
			fprintf(m_fp, "%g %g %g\n", dr.x, dr.y, dr.z);
		}		   

		list<string>::const_iterator it;
		for (it = point_data_fields.begin(); it != point_data_fields.end(); ++it) {

			printf("writing %s\n",it->c_str());

			char szname[256];
			strcpy(szname, it->c_str());
			Space2_(szname);

			// create the plot variable
			FEPlotData* pd = fecore_new<FEPlotData>(it->c_str(), m_fem);
			printf("region type %d\n", pd->RegionType());
			assert(pd->RegionType() == FE_REGION_NODE);
			int ndata = pd->VarSize(pd->DataType());
			printf("%d\n", ndata);

			// loop over all node sets
            // right now there is only one, namely the node set of all mesh nodes
            // so we just pass the mesh
			FEDataStream val; val.reserve(ndata * N);
			if (pd->Save(m_fem->GetMesh(), val))
			{

				printf("val size=%d", (int)val.size());

				// write the value array
				int ntype = pd->DataType();
				if (ntype == PLT_FLOAT) {
					fprintf(m_fp, "%s %s %s\n", "SCALARS", szname, "float");
					fprintf(m_fp, "%s %s\n", "LOOKUP_TABLE", "default");
					for (int i = 0; i<val.size(); ++i) fprintf(m_fp, "%g\n", val[i]);
				}
				else if (ntype == PLT_VEC3F) {
					fprintf(m_fp, "%s %s %s\n", "VECTORS", szname, "float");
					for (int i = 0; i<val.size(); i += 3) fprintf(m_fp, "%g %g %g\n", val[i], val[i + 1], val[i + 2]);
				}
				else if (ntype == PLT_MAT3FS) {
					fprintf(m_fp, "%s %s %s\n", "TENSORS", szname, "float");
					for (int i = 0; i<val.size(); i += 6)
						fprintf(m_fp, "%g %g %g\n%g %g %g\n%g %g %g\n\n",
							val[i    ], val[i + 3], val[i + 5],
							val[i + 3], val[i + 1], val[i + 4],
							val[i + 5], val[i + 4], val[i + 2]);
				}
				else if (ntype == PLT_MAT3FD) {
					fprintf(m_fp, "%s %s %s\n", "TENSORS", szname, "float");
					for (int i = 0; i<val.size(); i += 3)
						fprintf(m_fp, "%g %g %g\n%g %g %g\n%g %g %g\n\n",
							val[i], 0.f, 0.f,
							0.f, val[i + 1], 0.f,
							0.f, 0.f, val[i + 2]);
				}
			}
		}
	//}
}

//-----------------------------------------------------------------------------
void FEVTKExport::WriteCellData()
{

	FEMesh& mesh = m_fem->GetMesh();

	fprintf(m_fp, "\nCELL_DATA %d\n", mesh.Elements());
	list<string>::const_iterator it;
	for (it = cell_data_fields.begin(); it != cell_data_fields.end(); ++it) {

		char szname[256];
		strcpy(szname, it->c_str());
		Space2_(szname);

		// create the plot variable
		FEPlotData* pd = fecore_new<FEPlotData>(it->c_str(), m_fem);
		int ndata = pd->VarSize(pd->DataType());
		int ntype = pd->DataType();
		printf("region type %d\n", pd->RegionType());
		assert(pd->RegionType() == FE_REGION_DOMAIN);
		printf("data type %d\n", ntype);

		// loop over all domains
		int NDOMS = mesh.Domains();
		for (int i = 0; i < NDOMS; ++i)
		{
			bool first = i == 0;  // write headers

			FEDomain& dom = mesh.Domain(i);

			int N = dom.Elements();

			FEDataStream val; val.reserve(ndata * N);

			if (pd->Save(dom, val))
			{
	 
				if (ntype == PLT_FLOAT) {
					if (first) {
						fprintf(m_fp, "%s %s %s\n", "SCALARS", szname, "float");
						fprintf(m_fp, "%s %s\n", "LOOKUP_TABLE", "default");
					}
					for (int i = 0; i < val.size(); ++i) fprintf(m_fp, "%g\n", val[i]);
				}
				else if (ntype == PLT_VEC3F) {
					if (first) {
						fprintf(m_fp, "%s %s %s\n", "VECTORS", szname, "float");
					}
					for (int i = 0; i < val.size(); i += 3) fprintf(m_fp, "%g %g %g\n", val[i], val[i + 1], val[i + 2]);
				}
				else if (ntype == PLT_MAT3FS) {
					if (first) {
						fprintf(m_fp, "%s %s %s\n", "TENSORS", szname, "float");
					}
					for (int i = 0; i < val.size(); i += 6)
						fprintf(m_fp, "%g %g %g\n%g %g %g\n%g %g %g\n\n",
							val[i], val[i + 3], val[i + 5],
							val[i + 3], val[i + 1], val[i + 4],
							val[i + 5], val[i + 4], val[i + 2]);
				}
				else if (ntype == PLT_MAT3FD) {
					if (first) {
						fprintf(m_fp, "%s %s %s\n", "TENSORS", szname, "float");
					}
					for (int i = 0; i < val.size(); i += 3)
						fprintf(m_fp, "%g %g %g\n%g %g %g\n%g %g %g\n\n",
							val[i], 0.f, 0.f,
							0.f, val[i + 1], 0.f,
							0.f, 0.f, val[i + 2]);
				}
			}
		}
	}
}