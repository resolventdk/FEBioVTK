#include "FEVTKExport.h"

#include <stdio.h>

#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEAnalysis.h>
#include "FECore/FECore.h"
#include "FECore/FEPlotData.h"
#include "FECore/FEDomain.h"
#include "FECore/log.h"

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkInformation.h>
#include <vtkSmartPointer.h>


FEVTKExport::FEVTKExport(FEModel* fem)
{
	m_fem = fem;
	//m_prefix = fem->GetTitle(); // 

	// initialize json based meta file following Paraview format
	// https://gitlab.kitware.com/paraview/paraview/blob/v5.5.0/Documentation/release/ParaView-5.5.0.md#json-based-new-meta-file-format-for-series-added
	
	FILE* fp = fopen((m_prefix + ".vtm.series").c_str(), "wt");
	fprintf(fp, "{\n");
	fprintf(fp, "  \"file-series-version\" : \"1.0\",\n");
	fprintf(fp, "  \"files\" : [\n");
	fclose(fp);

}


FEVTKExport::~FEVTKExport(void)
{
	// finalize json based meta file
	FILE* fp = fopen((m_prefix + ".vtm.series").c_str(), "a");
	fprintf(fp, "\n");  // finish last entry
	fprintf(fp, "  ]\n");
	fprintf(fp, "}\n");
	fclose(fp);

}

bool FEVTKExport::Save()
{

	// get current step and time
	int step = m_ndump;
	m_ndump++;  // increment number of dumped solutions
	double time = m_fem->GetTime().currentTime;

	// set vtk file name
	char szfile[255];
	sprintf(szfile, "%s%3.3d.vtm", m_prefix.c_str(), step);

	// add this step to the meta file
	FILE* fp = fopen((m_prefix + ".vtm.series").c_str(), "a");
	if (m_ndump > 1) fprintf(fp, ",\n");  // finish last entry
	fprintf(fp, "     { \"name\" : \"%s\", \"time\" : %g }", szfile, time);
	fclose(fp);

	// create MultiBlockDataSets

	// root
	vtkSmartPointer<vtkMultiBlockDataSet> root =
		vtkSmartPointer<vtkMultiBlockDataSet>::New();

	// domains
	vtkSmartPointer<vtkMultiBlockDataSet> domains =
		vtkSmartPointer<vtkMultiBlockDataSet>::New();
	root->SetBlock(0, domains);
	root->GetMetaData((unsigned int)0)->Set(vtkCompositeDataSet::NAME(), "Domains");

	// surfaces
	vtkSmartPointer<vtkMultiBlockDataSet> surfaces =
		vtkSmartPointer<vtkMultiBlockDataSet>::New();
	root->SetBlock(1, surfaces);
	root->GetMetaData((unsigned int)1)->Set(vtkCompositeDataSet::NAME(), "Surfaces");

	// add domains
	if (!AddDomains(domains)) return false;

	// add surfaces
	if (!AddSurfaces(surfaces)) return false;

	// write multiBlockDataSet using XML writer
	vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer =
	vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
	writer->SetFileName(szfile);
	writer->SetInputData(root);
	writer->Write();

	// all done!
	return true;
}


bool FEVTKExport::AddDomains(vtkSmartPointer<vtkMultiBlockDataSet> multiBlockDataSet)
{

	FEMesh& mesh = m_fem->GetMesh();

	// loop over all domains
	int NDOMS = mesh.Domains();
	for (int i = 0; i < NDOMS; ++i)
	{
		FEDomain& dom = mesh.Domain(i);

		// create an unstructured grid for this domain
		vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
			vtkSmartPointer<vtkUnstructuredGrid>::New();

		// convert domain nodes to vtk points and insert into vtk grid
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		for (int n = 0; n < dom.Nodes(); n++) {
			vec3d& r = dom.Node(n).m_r0;
			points->InsertNextPoint(r.x, r.y, r.z);
		}
		unstructuredGrid->SetPoints(points);

		// insert elements into vtk grid one by one
		vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
		for (int e = 0; e < dom.Elements(); ++e)
		{

			FEElement el = dom.ElementRef(e);

			// local (domain) connectivity
			vtkSmartPointer<vtkIdList> cell = vtkSmartPointer<vtkIdList>::New();
			vector<int>& enodes = el.m_lnode;
			for (std::vector<int>::iterator it = enodes.begin(); it != enodes.end(); ++it)
				cell->InsertNextId(*it);

			// insert as correct vtk type
			int vtk_type = -1;
			switch (el.Type()) {
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
			case FE_HEX8RI: vtk_type = VTK_HEXAHEDRON; break;
			case FE_HEX8G1: vtk_type = VTK_HEXAHEDRON; break;
			case FE_HEX8G8: vtk_type = VTK_HEXAHEDRON; break;
			case FE_TET4G1: vtk_type = VTK_TETRA; break;
			case FE_TET4G4: vtk_type = VTK_TETRA; break;
			default: feLogEx(m_fem, "Unknown/unsupported volume element type!"); return false;  break;
			}

			unstructuredGrid->InsertNextCell(vtk_type, cell);

		} // elements, e


		// --- N O D E   D A T A ---

		// write displacement - not a FEPlotData object>??
		{
			vtkSmartPointer<vtkFloatArray> vtkArray =
				vtkSmartPointer<vtkFloatArray>::New();
			vtkArray->SetName("displacement");
			vtkArray->SetNumberOfComponents(3);
			vtkArray->SetNumberOfTuples(dom.Nodes());
			for (int n = 0; n < dom.Nodes(); n++) {
				vec3d dr = dom.Node(n).m_rt - dom.Node(n).m_r0;
				vtkArray->SetTuple3(n, dr.x, dr.y, dr.z);
			}
			unstructuredGrid->GetPointData()->AddArray(vtkArray);
		}

		for (list<string>::const_iterator var_name = point_data_fields.begin();
			var_name != point_data_fields.end(); ++var_name) {

			// create the plot variable
			FEPlotData* pd = fecore_new<FEPlotData>(var_name->c_str(), m_fem);

			// check that it is a valid nodal plot variable
			assert(pd->RegionType() == FE_REGION_NODE); // must be evaluated on the entire mesh :(

			// verify one value per node
			assert(pd->StorageFormat() == FMT_NODE);

			// get number of components
			int ndata = pd->VarSize(pd->DataType());

			// allocate data buffer
			FEDataStream val;
			val.reserve(ndata * dom.Nodes());

			// evaluate on the all nodes.. TODO: fix this
			if (pd->Save(mesh, val))
			{ 
				// create vtk array
				vtkSmartPointer<vtkFloatArray> vtkArray =
					vtkSmartPointer<vtkFloatArray>::New();
				vtkArray->SetName(var_name->c_str());
				vtkArray->SetNumberOfComponents(ndata);

				// transfer from buffer
				switch (pd->DataType()) {
				case PLT_FLOAT:  // scalar
					for (int j = 0; j < dom.Nodes(); j++) {
						int k = dom.NodeIndex(j); // global index of local node
						vtkArray->InsertNextTuple1(val[k]);
					}
					break;
				case PLT_VEC3F:  // 3 comp vector
					vtkArray->SetComponentName(0, "X");
					vtkArray->SetComponentName(1, "Y");
					vtkArray->SetComponentName(2, "Z");
					for (int j = 0; j < dom.Nodes(); j++) {
						int k = 3 * dom.NodeIndex(j); // global index of local node
						vtkArray->InsertNextTuple3(
							val[k], val[k + 1], val[k + 2]
						);
					}
					break;
				case PLT_MAT3FS:  // symmetric 3x3 matrix
					vtkArray->SetComponentName(0, "XX");
					vtkArray->SetComponentName(1, "YY");
					vtkArray->SetComponentName(2, "ZZ");
					vtkArray->SetComponentName(3, "XY");
					vtkArray->SetComponentName(4, "YZ");
					vtkArray->SetComponentName(5, "XZ");

					for (int j = 0; j < dom.Nodes(); j++) {
						int k = 6 * dom.NodeIndex(j); // global index of local node
						vtkArray->InsertNextTuple6(
							val[k], val[k + 1], val[k + 2],
							val[k + 3], val[k + 4], val[k + 5]
						);
					}
					break;
				case PLT_MAT3FD: // diagonal 3x3 matrix
					vtkArray->SetComponentName(0, "XX");
					vtkArray->SetComponentName(1, "YY");
					vtkArray->SetComponentName(2, "ZZ");
					for (int j = 0; j < dom.Nodes(); j++) {
						int k = 3 * dom.NodeIndex(j); // global index of local node
						vtkArray->InsertNextTuple3(
							val[k], val[k + 1], val[k + 2]
						);
					}
					break;
				default:
					feLogEx(m_fem, "Unknown/unsupported data type of plot variable!");
					return false;
					break;
				}

				// add to grid
				unstructuredGrid->GetPointData()->AddArray(vtkArray);

			} else {
				feLogEx(m_fem, "Unable to write node data: '%s'\n", var_name->c_str());
			}

		}  // var_name iterator


		// --- E L E M E N T   C E L L   D A T A ---

		// Write Material ID, then we can identify rigid domains
		{
			vtkSmartPointer<vtkIntArray> vtkArray =
				vtkSmartPointer<vtkIntArray>::New();
			vtkArray->SetName("MatID");
			vtkArray->SetNumberOfComponents(1);
			
			for (int e = 0; e < dom.Elements(); ++e)
			{
				FEElement el = dom.ElementRef(e);
				vtkArray->InsertNextTuple1(el.GetMatID());
			}
			
			unstructuredGrid->GetCellData()->AddArray(vtkArray);
		}

		// data fields
		for (list<string>::const_iterator var_name = cell_data_fields.begin();
			var_name != cell_data_fields.end(); ++var_name) {

			// create the plot variable
			FEPlotData* pd = fecore_new<FEPlotData>(var_name->c_str(), m_fem);

			// verify domain plot variable
			assert(pd->RegionType() == FE_REGION_DOMAIN);

			// verify one value per element
			assert(pd->StorageFormat() == FMT_ITEM);

			// get number of components
			int ndata = pd->VarSize(pd->DataType());

			// allocate data buffer
			FEDataStream val;
			val.reserve(ndata * dom.Elements());

			// evaluate on the mesh nodes
			// may fail if material does not support plotdata, eg. rigid material does not compute stresses
			if (pd->Save(dom, val))
			{
				// create vtk array
				vtkSmartPointer<vtkFloatArray> vtkArray =
					vtkSmartPointer<vtkFloatArray>::New();
				vtkArray->SetName(var_name->c_str());
				vtkArray->SetNumberOfComponents(ndata);

				// transfer from buffer		
				switch (pd->DataType()) {
				case PLT_FLOAT:  // scalar
					for (int j = 0; j < val.size(); j += 1)
						vtkArray->InsertNextTuple1(val[j]);
					break;
				case PLT_VEC3F:  // 3 comp vector
					vtkArray->SetComponentName(0, "X");
					vtkArray->SetComponentName(1, "Y");
					vtkArray->SetComponentName(2, "Z");
					for (int j = 0; j < val.size(); j += 3)
						vtkArray->InsertNextTuple3(
							val[j], val[j + 1], val[j + 2]
						);
					break;
				case PLT_MAT3FS:  // symmetric 3x3 matrix
					vtkArray->SetComponentName(0, "XX");
					vtkArray->SetComponentName(1, "YY");
					vtkArray->SetComponentName(2, "ZZ");
					vtkArray->SetComponentName(3, "XY");
					vtkArray->SetComponentName(4, "YZ");
					vtkArray->SetComponentName(5, "XZ");
					for (int j = 0; j < val.size(); j += 6)
						vtkArray->InsertNextTuple6(
							val[j], val[j + 1], val[j + 2],
							val[j + 3], val[j + 4], val[j + 5]
						);
					break;
				case PLT_MAT3FD: // diagonal 3x3 matrix
					vtkArray->SetComponentName(0, "XX");
					vtkArray->SetComponentName(1, "YY");
					vtkArray->SetComponentName(2, "ZZ");
					for (int j = 0; j < val.size(); j += 3)
						vtkArray->InsertNextTuple3(
							val[j], val[j + 1], val[j + 2]
						);
					break;
				default:
					feLogEx(m_fem, "Unknown/unsupported data type of plot variable!");
					return false;
					break;
				}

				// add to grid
				unstructuredGrid->GetCellData()->AddArray(vtkArray);
			}
			else {
				feLogEx(m_fem, "Unable to write element data '%s' for domain '%s'\n", var_name->c_str(), dom.GetName().c_str());
			}

		} // var_name iterator

		multiBlockDataSet->SetBlock(i, unstructuredGrid);
		multiBlockDataSet->GetMetaData((unsigned int)i)->Set(vtkCompositeDataSet::NAME(), dom.GetName().c_str());

	} // domain, i

	return true;

}


bool FEVTKExport::AddSurfaces(vtkSmartPointer<vtkMultiBlockDataSet> multiBlockDataSet)
{

	// Each surface is stored in its own file.
	// Surfaces are duplicated whenever referenced in multiple contact-pairs.
	// Merge dublicate surfaces before export.

    // get "surfaces"
	FEMesh& m = m_fem->GetMesh();
	int NS = m.Surfaces();

	// first create a list of surface names 
	// because surfaces are duplicated whenever referenced in a surface?/contact?-pair  
	std::vector<std::string> names;
	for (int i = 0; i < NS; ++i)
		names.push_back(m.Surface(i).GetName());

	// find unique names
	std::vector<std::string> unique_names;
	std::sort(names.begin(), names.end());  // need to be sorted before using unique
	std::unique_copy(names.begin(), names.end(), std::back_inserter(unique_names));

	// loop over unique surface names
	for (int i = 0; i < unique_names.size(); ++i)
	{

		// get the name
		const std::string name = unique_names[i];

		// get reference surface, count elements and size of data
		FESurface* Si = m.FindSurface(name); // what surface will this return???

		// write points
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		for (int n = 0; n < Si->Nodes(); n++)
		{
			FENode node = Si->Node(n);
			points->InsertNextPoint(node.m_r0.x, node.m_r0.y, node.m_r0.z);
		}

		// create a polydata object and add points to it
		vtkSmartPointer<vtkPolyData> polydata =
			vtkSmartPointer<vtkPolyData>::New();
		polydata->SetPoints(points);

		// create elements one by one
		vtkSmartPointer<vtkCellArray> cells =
			vtkSmartPointer<vtkCellArray>::New();

		for (int e = 0; e < Si->Elements(); e++)
		{
			FEElement el = Si->Element(e);
			vector<int>& enodes = el.m_lnode;

			vtkSmartPointer<vtkIdList> cell = vtkSmartPointer<vtkIdList>::New();
			for (std::vector<int>::iterator it = enodes.begin(); it != enodes.end(); ++it)
				cell->InsertNextId(*it);

			// insert as correct vtk type
			int vtk_type = -1;
			switch (el.Type()) {
				//FE_QUAD4G4,
				//FE_QUAD4NI,
				//FE_TRI3G1,
				//FE_TRI3G3,
				//FE_TRI3G7,
				//FE_TRI3NI,
				//FE_TRI6G3,
				//FE_TRI6G4,
				//FE_TRI6G7,
				////	FE_TRI6MG7,
				//FE_TRI6GL7,
				//FE_TRI6NI,
				//FE_TRI7G3,
				//FE_TRI7G4,
				//FE_TRI7G7,
				//FE_TRI7GL7,
				//FE_TRI10G7,
				//FE_TRI10G12,
				//FE_QUAD8G9,
				//FE_QUAD8NI,
				//FE_QUAD9G9,
				//FE_QUAD9NI,
			case FE_QUAD4G4: vtk_type = VTK_QUAD; break;
			case FE_QUAD4NI: vtk_type = VTK_QUAD; break;
			case FE_TRI3G1: vtk_type = VTK_TRIANGLE; break;
			case FE_TRI3G3: vtk_type = VTK_TRIANGLE; break;
			case FE_TRI3G7: vtk_type = VTK_TRIANGLE; break;
			case FE_TRI3NI: vtk_type = VTK_TRIANGLE; break;
			default: feLogEx(m_fem, "Unknown/unsupported surface element type!"); return false;  break;
			}

			// inserting directly into polydata cause seg error (ulike for vtugrid)
			cells->InsertNextCell(cell);  
		}
		polydata->SetPolys(cells);

		// --- S U R F A C E   P O I N T   D A T A ---
		{
			vtkSmartPointer<vtkFloatArray> vtkArray =
				vtkSmartPointer<vtkFloatArray>::New();
			vtkArray->SetName("displacement");
			vtkArray->SetNumberOfComponents(3);
			vtkArray->SetNumberOfTuples(Si->Nodes());

			for (int n = 0; n < Si->Nodes(); n++)
			{
				FENode node = Si->Node(n);
				vec3d dr = node.m_rt - node.m_r0;
				vtkArray->SetTuple3(n, dr.x, dr.y, dr.z);
			}
			polydata->GetPointData()->AddArray(vtkArray);
		}

		// --- S U R F A C E   C E L L   D A T A ---
		{
			vtkSmartPointer<vtkFloatArray> vtkArray =
				vtkSmartPointer<vtkFloatArray>::New();
			vtkArray->SetName("surface normal");
			vtkArray->SetNumberOfComponents(3);
			vtkArray->SetNumberOfTuples(Si->Elements());

			// use normal in displaced element (natural coord (0.5,0.5) is center)					
			for (int e = 0; e < Si->Elements(); e++)
			{
				vec3d n = Si->SurfaceNormal(Si->Element(e), 0.5, 0.5);
				vtkArray->SetTuple3(e, n.x, n.y, n.z);
			}
			polydata->GetCellData()->AddArray(vtkArray);
		}

		for (list<string>::const_iterator var_name = surface_data_fields.begin();
			var_name != surface_data_fields.end(); ++var_name) {

			// create the plot variable
			FEPlotData* pd = fecore_new<FEPlotData>(var_name->c_str(), m_fem);

			// verify region
			assert(pd->RegionType() == FE_REGION_SURFACE);

			// verify..? cell data not point data?
			assert(pd->StorageFormat() == FMT_ITEM);

			// get number of components
			int ndata = pd->VarSize(pd->DataType());

			// total size of data
			int nsize = ndata * Si->Elements();

			// allocate buffer and initialize to zero
			double* merge_val = new double[nsize];
			for (int k = 0; k < nsize; ++k) {
				merge_val[k] = 0.0;
			}

			// loop over all surfaces, match name and merge
			for (int j = 0; j < NS; ++j)
			{

				FESurface& Sj = m.Surface(j);

				// continue if names does not match
				if (name.compare(Sj.GetName()) != 0)
					continue;

				// check that number of elements matches
				assert(Si->Elements() == Sj.Elements());

				// get data for this surface
				FEDataStream val; val.reserve(nsize);
				if (pd->Save(Sj, val))
				{
					// merge using addition
                    // TODO: for some variables ie contact gap we should take min value
					for (int k = 0; k < nsize; ++k) {
						merge_val[k] += val[k];
					}
				}
				else {
					feLogEx(m_fem, "Unable to write surface data '%s' for surface '%s'\n", var_name->c_str(), name.c_str());
				}

			} // Surface(j), j

			// create vtk array
			vtkSmartPointer<vtkFloatArray> vtkArray =
				vtkSmartPointer<vtkFloatArray>::New();
			vtkArray->SetName(var_name->c_str());
			vtkArray->SetNumberOfComponents(ndata);

			// transfer from buffer
			switch (pd->DataType()) {
			case PLT_FLOAT:  // scalar
				for (int j = 0; j < nsize; j += 1)
					vtkArray->InsertNextTuple1(merge_val[j]);
				break;
			case PLT_VEC3F:  // 3 comp vector
				vtkArray->SetComponentName(0, "X");
				vtkArray->SetComponentName(1, "Y");
				vtkArray->SetComponentName(2, "Z");
				for (int j = 0; j < nsize; j += 3)
					vtkArray->InsertNextTuple3(
						merge_val[j], merge_val[j + 1], merge_val[j + 2]
					);
				break;
			default:
				feLogEx(m_fem, "Unknown/unsupported data type of plot variable!");
				return false;
				break;
			}

			// add to grid
			polydata->GetCellData()->AddArray(vtkArray);

			// clean up buffer
			delete merge_val;

		}  // plot variable

		multiBlockDataSet->SetBlock(i, polydata);
		multiBlockDataSet->GetMetaData((unsigned int)i)->Set(vtkCompositeDataSet::NAME(), name.c_str());

	} // unique surface names, i

	return true;

}
