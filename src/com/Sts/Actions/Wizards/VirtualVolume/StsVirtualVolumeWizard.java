package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Types.StsSeismicBoundingBox;

public class StsVirtualVolumeWizard extends StsWizard
{
    private StsMathVVolume mathVirtualVolume;
    private StsBlendVVolume blendVirtualVolume;
    private StsCrossplotVVolume xplotVirtualVolume;
    private StsRgbaVVolume rgbaVirtualVolume;
    private StsFilterVVolume filterVirtualVolume;
    private StsEPFVVolume epfVirtualVolume;
    private StsTypeVirtualVolume typeVirtualVolume;
    private StsDefineFloatVolume defineFloatVolume;
    private StsVirtualVolumeProcess processVolume;

    String name = "VVName";
    boolean floatVolume = false;
    public Object[] selectedSensors = new Object[0];
    
    private StsWizardStep[] mySteps =
    {
        typeVirtualVolume = new StsTypeVirtualVolume(this),
        mathVirtualVolume = new StsMathVVolume(this),
        blendVirtualVolume = new StsBlendVVolume(this),
        xplotVirtualVolume = new StsCrossplotVVolume(this),
        rgbaVirtualVolume = new StsRgbaVVolume(this),
        filterVirtualVolume = new StsFilterVVolume(this),
        epfVirtualVolume = new StsEPFVVolume(this),
        defineFloatVolume = new StsDefineFloatVolume(this),
        processVolume = new StsVirtualVolumeProcess(this)
    };

    public StsVirtualVolumeWizard(StsActionManager actionManager)
    {
        super(actionManager, 400, 400);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Define Virtual Volume");
        if(!super.start()) return false;

        checkAddToolbar();

        return true;
    }

    public boolean end()
    {
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
    	if(currentStep == typeVirtualVolume)
    	{
    		gotoTypeStep();
    	}
    	else if(currentStep == defineFloatVolume)
    		gotoStep(processVolume);
    	else
    	{
    		if(isFloatVolume())
    			gotoStep(defineFloatVolume);
    		else
    			gotoStep(processVolume);
    	}
    }
    
    public boolean isFloatVolume()
    {
    	return floatVolume;
    }
    
    public void setFloatVolume(boolean isFloat)
    {
    	floatVolume = isFloat;   	
    }
    
    public void gotoTypeStep()
    {
        switch(typeVirtualVolume.getVolumeType())
        {
            case StsVirtualVolume.SEISMIC_MATH:
                gotoStep(mathVirtualVolume);
                break;
            case StsVirtualVolume.SEISMIC_BLEND:
                gotoStep(blendVirtualVolume);
                break;
            case StsVirtualVolume.SEISMIC_XPLOT_MATH:
                if(getCrossplots().length == 0)
                    new StsMessage(this.frame,StsMessage.WARNING,"No crossplots currently available in project.");
                else
                    gotoStep(xplotVirtualVolume);
                break;
            case StsVirtualVolume.RGB_BLEND:
                gotoStep(rgbaVirtualVolume);
                break;
            case StsVirtualVolume.SEISMIC_FILTER:
//                filterVirtualVolume.setVolumeName(typeVirtualVolume.getVolumeName());
                filterVirtualVolume.panel.setSeismicVolumes(getSeismicVolumes());
                filterVirtualVolume.panel.setVirtualVolumes(getVirtualVolumes());
                filterVirtualVolume.setVolumeName(typeVirtualVolume.getVolumeName());
                gotoStep(filterVirtualVolume);
                break;
            case StsVirtualVolume.EP_SEISMIC_FILTER:
            	epfVirtualVolume.panel.setSeismicVolumes(getSeismicVolumes());
            	epfVirtualVolume.panel.setVirtualVolumes(getVirtualVolumes());
            	epfVirtualVolume.setVolumeName(typeVirtualVolume.getVolumeName());
              gotoStep(epfVirtualVolume);
              break;
            default:
                break;
        }
    }
    
    public int getVolumeType()
    {
        return typeVirtualVolume.getVolumeType();
    }
    
    public void setVolumeName(String name) { this.name = name; }
    public String getVolumeName() { return name; }

    public StsSeismicVolume[] getSeismicVolumes()
    {
        return (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
    }
    public StsVirtualVolume[] getVirtualVolumes()
    {
        StsVirtualVolumeClass vvClass = (StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class);
        return vvClass.getVirtualVolumes();
    }
    public StsCrossplot[] getCrossplots()
    {
        return (StsCrossplot[])model.getCastObjectList(StsCrossplot.class);
    }

    public boolean buildVolume(StsProgressPanel panel)
    {	
    	if(isFloatVolume())
    	{
    		StsSeismicVolume[] volumes = getSeismicVolumes();
    		for(int i=0; i<volumes.length; i++)
    		{
    			if(volumes[i] instanceof StsVirtualVolume)
    			{
    				panel.appendLine("Cannot create float volume from virtual volume: " + volumes[i].getName());
    				panel.appendLine("Will create byte volume.");
    				floatVolume = false;
    				break;
    			}
    		}     
    	}
    	
        switch(typeVirtualVolume.getVolumeType())
        {
            case StsVirtualVolume.SEISMIC_MATH:
                mathVirtualVolume.buildVolume(panel);
                break;
            case StsVirtualVolume.SEISMIC_BLEND:
                blendVirtualVolume.buildVolume(panel);
                break;
            case StsVirtualVolume.SEISMIC_XPLOT_MATH:
                xplotVirtualVolume.buildVolume(panel);
                break;
            case StsVirtualVolume.RGB_BLEND:
            	//rgbaVirtualVolume.buildVolume(panel);
                break;
            case StsVirtualVolume.SEISMIC_FILTER:
                filterVirtualVolume.buildVolume(panel);
                break;
            case StsVirtualVolume.EP_SEISMIC_FILTER:
                epfVirtualVolume.buildVolume(panel);
                break;
            default:
                break;
        }
        return true;
    }
    
    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsVirtualVolumeWizard vvWizard = new StsVirtualVolumeWizard(actionManager);
        vvWizard.start();
    }
}
