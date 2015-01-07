package com.Sts.Actions.Wizards.AnalogueCube;

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

public class StsSelectVolume extends StsWizardStep
{
    StsSelectVolumePanel panel;
    StsHeaderPanel header;

    public StsSelectVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(450, 400));
        header.setTitle("Aspect Energy Analogue Cube Analysis");
        header.setSubtitle("Select Target PostStack3d to Analyze");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#AnalogueCube");
        header.setLogo("AspectLogo.gif");
        header.setInfoText(wizardDialog,"(1) Select the volume to use as target volume.\n" +
                           "(2) Constrain the target volume by selecting a subvolume to limit the analysis.\n" +
                           "(3) Specify a name for the resulting Analogue Cube.\n" +
                           "(4) Select whether the analysis is done on 8-bit or 32-bit datasets.\n" +
                           "  **** The default is 32-bit if the high resolution option was selected when\n" +
                           "       the SegY file was processed.\n" +
                           "(5) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/marketing/aspect/AspectAnalogueCube.html");
    }

    public StsSeismicVolume getVolume()
    {
        return panel.getVolume();
    }

    public StsSubVolume getSubVolume()
    {
        return panel.getSubVolume();
    }

    public String getOptionalName()
    {
        return panel.getOptionalName();
    }

    public boolean isDataFloat()
    {
        return panel.isDataFloat();
    }
    public boolean start()
    {
        panel.initialize();
		StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
		subVolumeClass.setIsVisible(false);
		subVolumeClass.setIsApplied(false);   // Does not appear to work for sssv
        return true;
    }

    public boolean end()
    {
        return true;
    }
}
