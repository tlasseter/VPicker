package com.Sts.Actions.Wizards.Seismic;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.*;

public class StsSeismicRowColRangeEdit extends StsWizardStep //implements Runnable
{
    StsSeismicRowColRangeEditPanel panel;
    StsHeaderPanel header;
    StsSeismicBoundingBox[] volumes;

    public StsSeismicRowColRangeEdit(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSeismicRowColRangeEditPanel("Range Edit Parameters");
        header = new StsHeaderPanel();
        setPanels(panel, header);

        header.setTitle("Range Edit");
        header.setSubtitle("Edit Output Volume Range");
        header.setInfoText(wizardDialog,"(1) Set the inline and crossline ranges and increment\n" +
                           "(2) Once all crop and data range values have been specified press the Next>> Button");
    }

    public void initialize(StsSeismicBoundingBox[] volumes)
    {
        panel.initialize(volumes);
    }

    public void initialize(StsSeismicBoundingBox volumes)
    {
        panel.initialize(new StsSeismicBoundingBox[] {volumes });
    }

    public void initializeExtend(StsSeismicBoundingBox volumes)
    {
        panel.initializeExtend(new StsSeismicBoundingBox[] {volumes });
    }

    public boolean start()
    {
        wizard.dialog.setTitle("Post-Stack Range Edit");
        return true;
    }

    public StsEditedBoundingBox getEditedBox()
    {
        return panel.getEditedBox();
    }
    public StsEditedBoundingBox[] getCroppedBoxes()
    {
        return panel.getEditedBoxes();
    }

    public boolean end()
    {
        boolean success = true;
        return success;
    }
}