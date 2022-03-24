#include <FECore/sdk.h>
#include <FECore/Callback.h>
#include <FECore/FEModel.h>
#include <FECore/FECallBack.h>
#include "FECore/log.h"

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
        FEModel* pfem = &fem;

        if (nwhen == CB_SOLVED)
        {        
            feLogEx(pfem, "Writing VTK file\n");
            feLogEx(pfem, "===========================================================================\n");
            if (!vtkExporter->Save())
                feLogEx(pfem, "Failed to write VTK file!\n");
            feLogEx(pfem, "Done\n\n");

            delete vtkExporter; // plugin deconstructor not called on normal termination! 
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
