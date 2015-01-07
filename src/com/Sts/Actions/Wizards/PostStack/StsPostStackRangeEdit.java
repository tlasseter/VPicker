package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

public class StsPostStackRangeEdit extends StsWizardStep //implements Runnable
{
    StsPostStackRangeEditPanel panel;
    StsHeaderPanel header;
    boolean cropped = true;

    public StsPostStackRangeEdit(StsWizard wizard)
    {
        super(wizard);
        panel = new StsPostStackRangeEditPanel((StsSeismicWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);

        header.setTitle("Range Edit");
        header.setSubtitle("Edit Output Volume Range");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ProcessSeismic");
        header.setInfoText(wizardDialog,"(1) Set the inline, crossline and time/depth limits to crop\n" +
                           "(2) Set the data scaling limits. The maximum and minumum supplied were generated" +
                           " via a random scan small percentage of the data.\n" +
                           "     ***** All data will be converted to 8bit/32bit after scaled to user specified limits *****\n" +
                           "(3) Decide whether to set clipped values to the maxima or to null by setting the Null Clip checkbox\n" +
                           "(4) Once all crop and data range values have been specified press the Next>> Button");
    }

    public boolean start()
    {
        wizard.dialog.setTitle("Post-Stack Range Edit");
        try
        {
            panel.initialize();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "start", e);
            return false;
        }
        return true;
    }

    public boolean end()
    {
        boolean success = true;
        if(cropped && ((StsSeismicWizard) wizard).cropEnabled())
        {
            StsSeismicBoundingBox[] volumes = ((StsSeismicWizard) wizard).getSegyVolumesToProcess();
            StsCroppedBoundingBox croppedBox = panel.cropped3dRangeBox.croppedBoundingBox;
            for(int i=0; i<volumes.length; i++)
            {
                StsCroppedBoundingBox croppedBoxCopy = new StsCroppedBoundingBox(volumes[i], croppedBox, false);
                ((StsSegyVolume)volumes[i]).setCropBox(croppedBoxCopy);
            }
        }
        return success;
    }
}
