package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineWellSubVolume extends StsWizardStep
{
    StsDefineWellSubVolumePanel panel;
    StsHeaderPanel header;

    public StsDefineWellSubVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineWellSubVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Sub-Volume Definition");
        header.setSubtitle("Defining Well Constrained Sub-Volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SubVolume");
        header.setInfoText(wizardDialog,"(1) Specify the name of the well set sub-volume.\n" +
                           "(2) Specify the default cylinder color.\n" +
                           "(3) Specify the default cylinder height and radius.\n" +
                           "(4a) Select the Create Button to enter create mode and start creating well cylinders.\n" +
                           "    **** Graphically select the well to create the well cylinder ****\n" +
                           "(4b) Select the Edit Button to enter edit mode and edit existing well cylinders.\n" +
                           "    **** Graphically select the top or bottom handle on the cylinder to re-size it.\n" +
                           "    **** Graphically move the cylinder up or down the well path.\n" +
                           "(4c) Select the Delete Button to enter delete mode and delete existing cylinders.\n" +
                           "    **** Graphically select the center of the cylinder to delete ****\n" +
                           "(5) Press the Next> Button to move to the sub-volume setup screen.\n");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public StsWell getWell()
    {
        return panel.getWell();
    }
    public float getTopZ()
    {
        return panel.getTopZ();
    }
    public float getBottomZ()
    {
        return panel.getBottomZ();
    }
    public float getRadius()
    {
        return panel.getRadius();
    }

    public boolean end()
    {
		model.subVolumeChanged();
        return true;
    }

    public String getSubVolumeName()
    {
        return panel.getSubVolumeName();
    }

}
