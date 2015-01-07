
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Horizons;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.StsProject;
import com.Sts.UI.*;

public class StsSelectSurfaces extends StsWizardStep
{
    StsOrderedListPanel panel;
    StsHeaderPanel header;
    StsSurface[] surfaces = null;
    int[] selected = null;

    public StsSelectSurfaces(StsWizard wizard)
    {
        super(wizard, new StsOrderedListPanel(), null, new StsHeaderPanel());
        panel = (StsOrderedListPanel) getContainer();
//        infoPanel = (JTextField) getInfoContainer();
//        infoPanel.setIcon(createIcon("selectSurfaces.gif"));
//		infoPanel.setText("Select surfaces\nto define\nhorizons.");
//        infoPanel.setPreferredSize(new Dimension(75,300));
//        panel.setTitle("Select surfaces to define horizons:");
    	panel.setListRenderer(new StsColorListRenderer());
//        panel.setPreferredSize(new Dimension(400,300));

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Horizon Construction");
        header.setSubtitle("Select Surfaces");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Horizons");
        header.setInfoText(wizardDialog,"(1) Select all surfaces to convert to horizons.\n" +
                                  "   ***** Horizons are on same grid ***** \n" +
                                  "   ***** Surfaces are not necessarily on same grid ***** \n" +
                                  "(2) Press the Next>> Button.");
    }

    public boolean start()
    {
        if(model == null) return true;

        // get the horizon surfaces
    	if(surfaces == null)
        {
            /* Moved check to wizard so GUI does not come up
			surfaces = (StsSurface[])model.getCastObjectList(StsSurface.class);
			if(surfaces == null || surfaces.length == 0)
            {

                new StsMessage(model.win3d, StsMessage.WARNING, "No imported surfaces available.");
                wizard.cancel();
                return false;
//                return true;
            }
            */

			surfaces = (StsSurface[])model.getCastObjectList(StsSurface.class);
            if(surfaces == null) return false;
            surfaces = StsSurface.getZDomainSurfaces(surfaces, model.getProject().getZDomain());
            if(surfaces.length == 0)
            {
                surfaces = (StsSurface[])model.getCastObjectList(StsSurface.class);
                if(StsYesNoDialog.questionValue(wizard.frame, "No surfaces exist in " + model.getProject().getZDomainString() + ".\nDo you want to switch domains?"))
                {
                    model.getProject().setZDomain(surfaces[0].getZDomainOriginal());
                    surfaces = StsSurface.getZDomainSurfaces(surfaces, model.getProject().getZDomain());                    
                }
            }
            StsSurface[] orderedSurfaces = (StsSurface[])StsSurface.sortSurfaces(surfaces);
            int nSurfaces = surfaces.length;

            StsColorListItem[] items = new StsColorListItem[nSurfaces];
            int[] selected = new int[nSurfaces];
            for(int n = 0; n < nSurfaces; n++)
            {
                StsSurface surface = orderedSurfaces[n];
                String surfaceName = surface.getName();
                items[n] = new StsColorListItem(surface.getStsColor(), surfaceName);
                selected[n] = n;
            }
            panel.setItems(items);
            if(selected != null) panel.setSelectedIndices(selected);
        }
        return surfaces != null;
    }

    public String trimName(String fullString)
    {
    	int index = fullString.indexOf("(");
        if( index > 0 )
			return fullString.substring(0, index-1);
        return fullString;
    }

    public boolean end() { return getSelectedSurfaces() != null; }

    public StsSurface[] getSelectedSurfaces()
    {
        StsSurface[] selectedSurfaces;
		Object[] items = panel.getSelectedItems();
       	int nItems = items == null ? 0 : items.length;
        if(nItems == 0) return new StsModelSurface[0];

        selectedSurfaces = new StsSurface[nItems];
        for( int i=0; i<nItems; i++ )
        {
            StsColorListItem item = (StsColorListItem)items[i];
            String name = trimName(item.getName());
            selectedSurfaces[i] = (StsSurface)StsMainObject.getListObjectWithName(surfaces, name);
        }
    	return selectedSurfaces;
    }
}
