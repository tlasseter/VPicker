package com.Sts.Actions.Wizards.AnalogueCube;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
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

public class StsDefineScope extends StsWizardStep
{
    StsDefineScopePanel panel;
    StsHeaderPanel header;

//    StsSeismicVolume sourceVolume = null;
//    StsBoxSubVolume sourceSubVolume = null;
//    StsBoxSetSubVolume boxSetSubVolume = null;
//    static StsBoxSetSubVolumeClass boxSetClass = null;
    final StsColor stsColor = StsColor.CYAN;

    public StsDefineScope(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineScopePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(450, 400));
        header.setTitle("Aspect Energy Analogue Cube Analysis");
        header.setSubtitle("Define Source SubVolume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#AnalogueCube");
        header.setLogo("AspectLogo.gif");
        header.setInfoText(wizardDialog,"(1) Define the source sub-volume using the mouse or the manual values.\n" +
                           " *** The box will update in the graphics window as text values are changed ***\n" +
                           "(2) Press the Next>> Button when done\n");
        header.setLink("http://www.s2ssystems.com/marketing/aspect/AspectAnalogueCube.html");
    }

    public boolean start()
    {
//        StsBoxSetSubVolumeClass boxSetClass = (StsBoxSetSubVolumeClass)model.getStsClass(StsBoxSetSubVolume.class);
 //       boxSetClass.setIsVisible(true);

        StsSeismicVolume sourceVolume = ((StsAnalogueCubeWizard)wizard).getTargetVolume();
        panel.initialize(sourceVolume);
        return true;
    }

    public boolean end()
    {
        if(panel.getSourceBox() != null)
        {
            // In case the user changed the volume after the box was defined.
            panel.defineSourceVolume();
            return true;
        }
        new StsMessage(model.win3d, StsMessage.WARNING, "Source box has not been selected in window yet.");
        return false;
    }

    public StsSeismicVolume getSourceVolume() { return panel.sourceVolume; }

    public StsBoxSubVolume getSourceBox() { return panel.getSourceBox(); }

	public StsBoxSetSubVolume getSourceBoxSet() { return panel.getSourceBoxSet(); }
    public void clearBox()
    {
        if(panel == null) return;
        panel.deleteBoxes();
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        return panel.performMouseAction(mouse, (StsGLPanel3d)glPanel);
    }
}
