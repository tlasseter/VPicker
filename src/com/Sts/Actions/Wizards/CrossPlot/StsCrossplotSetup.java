
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.CrossPlot;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Toolbars.*;

public class StsCrossplotSetup extends StsAction
{
    private StsModel model = null;
    private StsCrossplotToolbar crossplotToolbar;

    public StsCrossplotSetup(StsActionManager actionManager)
    {
        super(actionManager);
        model = actionManager.getModel();
        crossplotToolbar = new StsCrossplotToolbar(model.win3d);

        StsCrossplot[] crossplots = (StsCrossplot[])model.getCastObjectList(StsCrossplot.class);
        int nCrossplots = crossplots.length;
        StsCrossplotClass crossplotClass = StsCrossplot.getCrossplotClass();

        for(int n = 0; n < nCrossplots; n++)
        {
//            crossplotClass.setCurrentCrossplot(crossplots[n]);
            crossplots[n].initialize(model);
        }
        StsCrossplot currentCrossplot = crossplotClass.getCurrentCrossplot();
        if(currentCrossplot != null) currentCrossplot.processPolygons();
    }

    public boolean start() { return false; }

    public void checkAddToolbar()
    {
        if(!model.win3d.hasToolbarNamed(StsCrossplotToolbar.NAME))
        {
 //           StsViewSelectToolbar tb = model.win3d.getViewSelectToolbar();
 //           tb.setButtonVisibility(tb.getComponentNamed("XPView"), true);
            model.win3d.addToolbar(crossplotToolbar);
        }
    }
}

