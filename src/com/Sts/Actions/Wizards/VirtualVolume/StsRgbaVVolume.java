package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.*;
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

public class StsRgbaVVolume extends StsWizardStep
{
    StsRgbaVVolumePanel panel;
    StsHeaderPanel header;
	StsVirtualVolume virtualVolume = null;

    final static int RGBA = 0;
    final static int HSIA = 1;

    public StsRgbaVVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsRgbaVVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Virtual PostStack3d Definition");
        header.setSubtitle("Defining Color Virtual PostStack3d");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VirtualVolume");
        header.setInfoText(wizardDialog,"(1) Select a volume for each of the color components.\n" +
                           "     **** If no transparency volume is selected, all values will be opaque. ****\n" +
                           "(2) Press the Finish> Button to create the volume.\n" +
                           "     ***** A virtual volume combobox will be added to the toolbar panel to control the viewing *****\n");
    }

    public boolean start()
	{
        panel.initialize();
        wizard.enableFinish();
        disableNext();
        return true;
    }

    public boolean end()
    {
		if(virtualVolume != null) return true;

        String name = null;

        StsSeismicBoundingBox[] volumes = panel.getSelectedVolumes();
        //
        // Verify that the volumes have the same bounding box
        //
        if(volumes.length < 3)
        {
            new StsMessage(wizard.frame, StsMessage.ERROR,"Failed to construct virtual volume. Need at least three volumes.\n");
            return false;
        }
        for(int i=1; i<volumes.length; i++)
        {
            if(!volumes[i-1].sameAs(volumes[i]))
            {
                new StsMessage(wizard.frame, StsMessage.ERROR,"Failed to construct virtual volume. Volumes must have identical bounding box.\n");
                return false;
            }
        }
        //
        // Construct a name if one is not supplied
        //
        name = ((StsVirtualVolumeWizard)wizard).getVolumeName();
        if(name.equals("VVName") || (name == null) || (name.length() <= 0))
        {
            if(volumes.length == 3)
                name = makeVirtualVolumeName(panel.getColorType(), volumes[0].getName(),
                                             volumes[1].getName(), volumes[2].getName(), null);
            else
                name = makeVirtualVolumeName(panel.getColorType(), volumes[0].getName(),
                                             volumes[1].getName(), volumes[2].getName(), volumes[3].getName());
            if(name == null)
            {
                new StsMessage(wizard.frame, StsMessage.ERROR,"Failed to construct virtual volume. Logical name invalid.\n");
                return false;
            }
        }
        virtualVolume = new StsRGBAVirtualVolume(volumes, name);
        if (virtualVolume == null)
        {
            new StsMessage(wizard.frame, StsMessage.ERROR,"Failed to construct virtual volume.\n");
            return false;
        }
    	virtualVolume.computeHistogram();        
        return true;
    }


    public String makeVirtualVolumeName(int type, String v1, String v2, String v3, String v4)
    {
       String operString = null;
       String l1,l2,l3,l4;

       if(type == RGBA)
           {l1 = "R_"; l2 = "G_"; l3 = "B_"; l4 = "A_";}
       else
           {l1 = "H_"; l2 = "S_"; l3 = "I_"; l4 = "A_";}

       if(v1.length() > 5)
           v1 = l1 + v1.substring(0,5);
       else
           v1 = l1 + v1;
       if(v2.length() > 5)
           v2 = l2 + v2.substring(0,5);
       else
           v2 = l2 + v2;
       if(v3.length() > 5)
           v3 = l3 + v3.substring(0,5);
       else
           v3 = l3 + v3;
       if(v4 != null)
       {
           if(v4.length() > 5)
               v4 = l4 + v4.substring(0, 5);
           else
               v4 = l4 + v4;
       }
       return v1 + v2 + v3 + v4;
    }
}
