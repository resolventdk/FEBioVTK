#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>

#include <FECore/sdk.h>
#include <FECore/Callback.h>
#include <FECore/FEModel.h>
#include <FECore/FECallBack.h>
#include "FEVTKExport.h"

class MyCallback : public FECallBack
{
private:
    FEVTKExport* vtkExporter;

public:
    MyCallback(FEModel* pfem) : FECallBack(pfem, CB_ALWAYS) 
    {
        vtkExporter = new FEVTKExport(pfem);
    }

    ~MyCallback()
    {
        delete vtkExporter;
    }

    bool Execute(FEModel& fem, int nwhen)
    {
		if (nwhen == CB_MAJOR_ITERS)
		{
            printf("*hjs: Writing VTK file\n");
            if (!vtkExporter->Save())
				printf("*hjs: Failed to write vtk file!\n");
			printf("*hjs: All done\n");

        }

        if (nwhen == CB_SOLVED)
        {
            delete vtkExporter;
        }

		// all done!
        return true;
    }

};

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);

	// feature classes
    REGISTER_FECORE_CLASS(MyCallback, "vtk_callback");

}

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}
