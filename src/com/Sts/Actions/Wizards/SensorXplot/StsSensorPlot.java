package com.Sts.Actions.Wizards.SensorXplot;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsSensorPlot extends StsWizardStep
{
	StsSensorPlotPanel panel;
    StsHeaderPanel header;

    public StsSensorPlot(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSensorPlotPanel((StsSensorXplotWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Crossplot Sensors");
        header.setSubtitle("Define Plot");
        header.setInfoText(wizardDialog,"(1) Select the desired attribute for the X-axis.\n" +
        		"(2) Select the desired attribute for the Y-axis\n" +
        		" **** Each time an attribute is changed the polygons are cleared.\n" +
        		"(3) Select the desired polygon color.\n" +
        		" **** One polygon is allowed per color.\n" +
        		"(4) Draw the polygon on the plot by clicking on the desired perimeter points.\n" +
        		"(5) Select whether to plot events inside or outside the drawn polygons\n" +
        		"(6) Press the Clear All or Clear Current to remove all polygons or the current polygon.\n" +
        		"(7) Add another level to the crossplot analysis\n" +
        		" **** Only the events inside polygons will be included in next level.\n" + 
        		"(8) Return to the sensor selection step if additional sensors are to be crossplotted\n" +
        		"(9) Current results can be exported by selecting the sensor folder with the middle button on the object panel.\n" +
        		"(10) Press the Finish Button, when done to clear the crossplot and reset the 3D view.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorXplot");
    }

    public boolean start()
    {
        panel.initialize();
        wizard.enableFinish();
        return true;
    }
    
    public boolean end()
    {
        return true;
    }

    public StsColor getClusterStsColor(int index)
    {
        return panel.getClusterStsColor(index);
    }
}
