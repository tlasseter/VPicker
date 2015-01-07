package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Progress.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsBlendVVolume extends StsWizardStep
{
    StsBlendVVolumePanel panel;
    StsHeaderPanel header;
	StsBlendedVirtualVolume virtualVolume = null;

    public StsBlendVVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsBlendVVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Virtual Volume Definition");
        header.setSubtitle("Defining Blended Virtual Volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VirtualVolume");
        header.setInfoText(wizardDialog,"(1) Select the volumes to blend.\n" +
                           "(2) Select the logical operation to apply.\n" +
                           "(3) Press the Finish> Button to create the volume.\n" +
                           "     ***** A virtual volume combobox will be added to the toolbar panel to control the viewing *****\n");

    }

    public boolean start()
	{
        panel.initialize();
        return true;
    }

    public boolean end()
    {
    	return true;
    }
    
    public boolean buildVolume(StsProgressPanel ppanel)
    {
 		if(virtualVolume != null) return true;
        String name = null;

        StsSeismicVolume[] volumes = panel.getSelectedVolumes();
        //
        // Verify that the volumes have the same bounding box
        //
        if (volumes.length == 2)
        {
            if (!volumes[0].sameAs(volumes[1]))
            {
                ppanel.appendLine("Failed to construct virtual volume. Volumes must have identical bounding box.\n");
                wizard.gotoStep(this);
            }
        }
        double condition = panel.getScalar();
        int logical = panel.getOperator();
        //
        // Construct a name if one is not supplied
        //
        name = ((StsVirtualVolumeWizard)wizard).getVolumeName();
        if(name.equals("VVName") || (name == null) || (name.length() <= 0))
        {
            if(volumes.length == 2)
                name = makeVirtualVolumeName(volumes[0].getName(),  volumes[1].getName(),
                                             new Double(condition).toString(), logical);
            else
                name = makeVirtualVolumeName(volumes[0].getName(), null, new Double(condition).toString(), logical);
            if(name == null)
            {
            	ppanel.appendLine("Failed to construct blended virtual volume. Logical name invalid.\n");
                wizard.gotoStep(this);
            }
        }
        ppanel.appendLine("Constructing blended virtual volume:" + name);
        virtualVolume = new StsBlendedVirtualVolume(volumes, name, condition, logical);
        if (virtualVolume == null)
        {
        	ppanel.appendLine("Failed to construct blended virtual volume.\n");
            return false;
        }
        if(((StsVirtualVolumeWizard)wizard).isFloatVolume())
        {
            ppanel.appendLine("Constructing float resolution blended volume:" + name);        	
    		StsVirtualVolumeConstructor.createVirtualVolume(model, volumes, virtualVolume, true, ppanel);
    	    success = true;
            ppanel.appendLine("Successfully created float resolution math volume:" + name);     	    
            wizard.enableFinish();
            model.win3dDisplay();
        } 
        else
        {
            virtualVolume.addToProject(false);
            virtualVolume.addToModel();
            ((StsBlendedVirtualVolume)virtualVolume).getBlendedVirtualVolumeClass().setIsVisibleOnCursor(true);
            virtualVolume.computeHistogram();
            model.viewObjectChangedAndRepaint(this, virtualVolume);
        }
        ppanel.appendLine("Successfully created blended virtual volume:" + name);
        return true;
    }

    public String makeVirtualVolumeName(String v1, String v2, String condition, int logical)
    {
        String operString = null;
        if(v1.length() > 8)
            v1 = v1.substring(0,8);
        if(v2 != null)
            if(v2.length() > 8)
                v2 = v2.substring(0,8);
        switch(logical)
        {
            case StsBlendedVirtualVolume.GT:
                if(v2 != null)
                    operString = v1 + "_U_" + v2 + "GT" + condition;
                else
                    operString = v1 + "GT" + condition;
                break;
            case StsBlendedVirtualVolume.LT:
                if(v2 != null)
                    operString = v1 + "_U_" + v2 + "LT" + condition;
                else
                    operString = v1 + "LT" + condition;
                break;
            case StsBlendedVirtualVolume.LE:
                if(v2 != null)
                    operString = v1 + "_U_" + v2 + "LE" + condition;
                else
                    operString = v1 + "LE" + condition;
                break;
            case StsBlendedVirtualVolume.GE:
                if(v2 != null)
                    operString = v1 + "_U_" + v2 + "GE" + condition;
                else
                    operString = v1 + "GE" + condition;
                break;
            default:
                operString = null;
                break;
        }
        return operString;
    }

}
