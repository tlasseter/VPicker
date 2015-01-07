package com.Sts.Actions.Wizards.SurfaceCurvature;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author unascribed
 * @version 1.1
 */

public class StsPatchPick extends StsWizardStep
{
    public StsPatchPickPanel panel;
    public StsHeaderPanel header;
    private StsVolumeCurvatureWizard patchVolumeWizard;

    public StsPatchPick(StsWizard wizard)
    {
        super(wizard);
        patchVolumeWizard = (StsVolumeCurvatureWizard)wizard;

        panel = new StsPatchPickPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 600));
        header.setTitle("Pick Volume Definition");
        header.setSubtitle("Define Volume Picking Parameters");
        header.setInfoText(wizardDialog, "(1) Select the Pick Type.\n" +
                "   **** Maximun, Minimum or Min+Max. ***** \n" +
                "   **** The Selected event types will be used to pick patches.****\n" +
                "(2) Specify window size.\n" +
                "   **** Window Size is the correlation window to be tracked, ****\n" +
                "   **** i.e. the window to extract possible correlative events. ****\n" +
                "(3) Specify Maximum difference between trace Picks.\n" +
                "   **** Pick Difference is effectively the allowable dip. ****\n" +
                "(4) Set the Minimum allowable patch size.\n" +
                "   **** Patches with fewer valid points than this will be eliminated. ****\n" +
                "(5) Set the Minimum allowable correlation value.\n" +
                "   **** Events that fail the correlation limit will not be added to the patch. ****\n" +
                "(6) Set the correlation wavelength multiplier.\n" +
                "   **** Trace to trace correlation will use a window determined as the   ****\n" +
                "   **** wavelength multiplier times the dynamicly determined wavelength. ****\n" +
                "(5) Press the Next>> Button when picking is completed.");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeCurvature");
    }

    public boolean start()
    {
        panel.initialize();
        disableNext();
        // disableFinish();
        return true;
    }

    public void setStatusMessage(String text)
    {
        panel.setStatusLabel(text);
    }

    public boolean end()
    {
        return true;
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        return ((StsVolumeCurvatureWizard)wizard).patchPick(mouse, glPanel);
    }
/*
    public boolean performMouseActionPicking(StsMouse mouse, StsGLPanel glPanel)
    {

        StsPickItem pickItem;

        try
        {
            StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
            int leftButtonState = mouse.getButtonState(StsMouse.LEFT);

            if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
            {
                StsPoint cursorZPoint = model.win3d.getCursor3d().getPointInCursorPlane(glPanel3d, StsCursor3d.ZDIR, mouse);
                return true;
            }

            if(leftButtonState != StsMouse.RELEASED) return true;

            // mouseButton RELEASED: process selection

            mouse.clearButtonState(StsMouse.LEFT);
            StsMethod patchPickMethod = new StsMethod(patchVolumeWizard, "pickPatch", glPanel3d);
            if(!StsJOGLPick.pick3d(glPanel, patchPickMethod, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_FIRST)) return true;
            pickItem = StsJOGLPick.pickItems[0];
		    int patchID = pickItem.names[0];
            StsPatchVolume patchVolume = patchVolumeWizard.patchVolume;
            patchVolume.addSelectedPatch(patchID);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeOnSurface.performMouseAction() failed.",
                    e, StsException.WARNING);
            return false;
        }
    }
*/
}
