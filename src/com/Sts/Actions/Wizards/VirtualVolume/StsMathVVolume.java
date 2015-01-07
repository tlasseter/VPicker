package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsMathVVolume extends StsWizardStep
{
    StsMathVVolumePanel panel;
    StsHeaderPanel header;
	StsVirtualVolume virtualVolume = null;
	StsSeismicVolume[] volumes = null;
	
    public StsMathVVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsMathVVolumePanel((StsVirtualVolumeWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Virtual Volume Definition");
        header.setSubtitle("Defining Math Virtual Volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VirtualVolume");
        header.setInfoText(wizardDialog,"(1) Select the volumes and/or scalar to operate on.\n" +
                           "(2) Select the operation to be applied.\n" +
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
        volumes = panel.getSelectedVolumes();
        //
        // Verify that the volumes have the same bounding box
        //
        if (volumes.length == 2)
        {
            if (!volumes[0].sameAs(volumes[1]))
            {
                ppanel.appendLine("Failed to construct virtual volume. Volumes must have identical bounding box.");
                wizard.gotoStep(this);
            }
        }
        int operator = panel.getOperator();
        double scalar = panel.getScalar();
        if (operator == StsMathVirtualVolume.TRANSPARENCY)
        {
            ppanel.appendLine("Transparency is not enabled yet.");
            wizard.gotoStep(this);
            return false;
        }
        if (operator == StsVirtualVolume.UNDEFINED)
        {
            ppanel.appendLine("Failed to construct virtual volume. Must specify operator.");
            wizard.gotoStep(this);
            return false;
        }

        //
        // Construct a name if one is not supplied
        //
        name = ((StsVirtualVolumeWizard)wizard).getVolumeName();
        if(name.equals("VVName") || (name == null) || (name.length() <= 0))
        {
            if (volumes.length == 2)
                name = makeVirtualVolumeName(volumes[0].getName(), volumes[1].getName(), operator);
            else
                name = makeVirtualVolumeName(volumes[0].getName(), new Double(scalar).toString(), operator);
        }
        ppanel.appendLine("Constructing math virtual volume:" + name);
        
        virtualVolume = new StsMathVirtualVolume(volumes, name, operator, scalar);
        if (virtualVolume == null)
        {
            new StsMessage(wizard.frame, StsMessage.ERROR, "Failed to construct virtual volume.\n");
            return false;
        }
        ppanel.appendLine("Successfully created math virtual volume:" + name); 
        
        // Create concrete volume and remove virtual volume.
        if(((StsVirtualVolumeWizard)wizard).isFloatVolume())
        {
            ppanel.appendLine("Constructing float resolution math volume:" + name);        	
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
        	((StsMathVirtualVolume)virtualVolume).getMathVirtualVolumeClass().setIsVisibleOnCursor(true);
        	virtualVolume.computeHistogram();        	
        	model.viewObjectChangedAndRepaint(this, virtualVolume);
        }
        ppanel.appendLine("Press Finish>> Button to complete virtual volume creation.");
        return true;
    }

    public String makeVirtualVolumeName(String v1, String v2, int op)
    {
       String operString = null;
       if(v1.length() > 8)
           v1 = v1.substring(0,8);
       if(v2.length() > 8)
           v2 = v2.substring(0,8);
       switch(op)
       {
           case StsMathVirtualVolume.ADDITION:
               operString = v1 + "_pls_" + v2;
               break;
           case StsMathVirtualVolume.SUBTRACTION:
               operString = v1 + "_min_" + v2;
               break;
           case StsMathVirtualVolume.MULTIPLY:
               operString = v1 + "_mul_" + v2;
               break;
           case StsMathVirtualVolume.DIVIDE:
               operString = v1 + "_div_" + v2;
               break;
           case StsMathVirtualVolume.AVERAGE:
               operString = v1 + "_avg_" + v2;
               break;
           case StsMathVirtualVolume.MAXIMUM:
               operString = v1 + "_max_" + v2;
               break;
           case StsMathVirtualVolume.MINIMUM:
               operString = v1 + "_min_" + v2;
               break;
           case StsMathVirtualVolume.ABSOLUTE:
               operString = "abs_" + v1;
               break;
           case StsMathVirtualVolume.POWER:
               operString = "pow_" + v1 + "_to_" + v2;
               break;
           case StsMathVirtualVolume.LOG:
               operString = "log_" + v1;
               break;
           case StsMathVirtualVolume.ASQRT:
               operString = "abssqrt_" + v1;
               break;
           case StsMathVirtualVolume.COS:
               operString = "cos_" + v1;
               break;
           case StsMathVirtualVolume.SIN:
               operString = "sin_" + v1;
               break;
           case StsMathVirtualVolume.SQUARE:
               operString = v1 + "_squared";
               break;
           case StsMathVirtualVolume.MODULUS:
               operString = "mod_" + v1;
               break;
            case StsMathVirtualVolume.INVERSE:
                operString = "inv_" + v1;
                break;
            case StsMathVirtualVolume.TRANSPARENCY:
               operString = v1 + "_transparency_" + v2;
               break;
           case StsMathVirtualVolume.UNDEFINED:
           default:
               operString = v1 + "_undef_" + v2;
       }
       return operString;
   }

}
