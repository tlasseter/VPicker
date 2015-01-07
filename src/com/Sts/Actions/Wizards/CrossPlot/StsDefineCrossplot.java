package com.Sts.Actions.Wizards.CrossPlot;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;
import java.util.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineCrossplot extends StsWizardStep
{
    StsDefineCrossplotPanel panel;
    StsHeaderPanel header;
	private StsCrossplot crossplot = null;

    public StsDefineCrossplot(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineCrossplotPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 300));
        header.setTitle("Cross Plot Definition");
        header.setSubtitle("Defining Cross Plot");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Crossplot");
    }

    public int[] getAxisVolumeIndices()
    {
        return panel.getAxisVolumeIndices();
    }
    public int getAttributeVolumeIndex()
    {
        return panel.getAttributeVolumeIndex();
    }
    public boolean start()
    {
        panel.initialize();
		crossplot = null;
        return true;
    }

    public boolean end()
    {
		if(crossplot != null) return true;
		return buildCrossplot();
    }

	public boolean buildCrossplot()
	{
        StsSeismicVolume[] selectedVolumes = panel.getAxisSeismicVolumes();
        StsSeismicVolume attributeVolume = panel.getAttributeSeismicVolume();
        int nVolumes = selectedVolumes.length;
        if(nVolumes < 2) return false;
        if(!selectedVolumes[0].sameAs(selectedVolumes[1]))
        {
            new StsMessage(wizard.frame, StsMessage.ERROR, "The selected volumes are not the same, unable to crossplot.");
            return false;
        }
        for(int n = 0; n < nVolumes; n++)
        {
            if(selectedVolumes[n] == null)
            {
                new StsMessage(wizard.frame, StsMessage.ERROR, "No seismic volume selected for axis number " + n);
                return false;
            }
        }
        if(attributeVolume != null)
        {
            if(!selectedVolumes[0].sameAs(attributeVolume))
            {
                new StsMessage(wizard.frame, StsMessage.ERROR, "The selected volumes are not the same, unable to crossplot.");
                return false;
            }
            selectedVolumes = (StsSeismicVolume[])StsMath.arrayAddElement(selectedVolumes, attributeVolume);
        }
        
        StsCrossplotToolbar crossplotToolbar = (StsCrossplotToolbar)model.win3d.getToolbarNamed(StsCrossplotToolbar.NAME);
        if(crossplotToolbar != null)
        {
        	if(selectedVolumes.length == 3)
        		crossplotToolbar.enableDensityToggle(true);
        	else
        		crossplotToolbar.enableDensityToggle(false);
        }
        
        // Process Name
        String name = panel.getName();
        if(name.length() == 0)
        {
            name = "XP-" + selectedVolumes[0].getName().substring(0,8) + selectedVolumes[1].getName().substring(0,8);
            new StsMessage(wizard.frame, StsMessage.ERROR, "No name supplied, using default name: " + name);
        }

        try
        {
            model.disableDisplay();
            crossplot = StsCrossplot.constructor(name, selectedVolumes);
            crossplot.setTypeLibrary(((StsCrossplotWizard)wizard).getSelectedLibrary());

            // Add polygons if file was selected
            if(panel.getPolygonFile() != null)
                buildPolygonObjects(crossplot);

            model.enableDisplay();
            wizard.enableFinish();
        }
        catch(Exception e)
        {
            new StsMessage(wizard.frame, StsMessage.ERROR, "Failed to construct crossplot.\n" + e.getMessage());
            return false;
        }
        crossplot.addToModel();
        
        // Can only run wizard from main window
 //       StsViewXP viewXP = (StsViewXP)model.win3d.glPanel3d.checkAddView(StsViewXP.class);
//		StsViewXP viewXP = (StsViewXP)model.setViewPreferred(StsViewXP.class, StsViewXP.viewClassnameXP);
        if(!model.win3d.constructViewPanel(StsViewXP.viewNameXP, StsViewXP.class)) return false;
        StsViewXP viewXP = (StsViewXP) model.win3d.getGlPanel3d().view;
        viewXP.crossplotChanged = true;

        // in main window, toggle the XP view button
        // model.win3d.selectToolbarItem(StsViewSelectToolbar.NAME, viewXP.getViewName(), true);
        return true;
    }

    private boolean buildPolygonObjects(StsCrossplot crossplot)
    {
        StringTokenizer st = null;
        String line = null;
        String nameType = "none";
        StsType type = null;
        Integer nVertices = new Integer(0);
        float vertex[] = new float[2];
        float totalAxis[][] = null;
        int typeIdx = 0;
        String defaultLibraryName = new String("crossplot");

        try
        {
            BufferedReader polygonFile = new BufferedReader(new FileReader(panel.getPolygonFile()));
            while (true)
            {
                line = polygonFile.readLine();
                if (line == null)
                    break;
                if (line.equals("<Polygon>"))
                    break;
            }
            StsTypeLibrary typeLib = crossplot.getTypeLibrary();
            StsSpectrum basicSpectrum = model.getSpectrum("Basic");
            totalAxis = crossplot.getTotalAxisRanges();
            while (true)
            {
                line = polygonFile.readLine();
                if(line == null) break;
                st = new StringTokenizer(line);
                String dead = st.nextToken();
                StsColor color = new StsColor(new Float(st.nextToken()).floatValue(),
                                                 new Float(st.nextToken()).floatValue(),
                                                 new Float(st.nextToken()).floatValue());
                if(!dead.equals(nameType))
                {
                    nameType = dead;
                    type = crossplot.getPolygonType(nameType, color);
                    typeIdx++;
                }
                // Create and load polygon
                StsXPolygon polygon = new StsXPolygon(type);
                crossplot.addPolygon(polygon);
                nVertices = new Integer(st.nextToken());
                for(int i=0; i<nVertices.intValue(); i++)
                {
                    line = polygonFile.readLine();
                    st = new StringTokenizer(line);
                    vertex[0] = new Float(st.nextToken()).floatValue();
                    vertex[1] = new Float(st.nextToken()).floatValue();
                    if(vertex[0] < totalAxis[0][0]) vertex[0] = totalAxis[0][0];
                    if(vertex[0] > totalAxis[0][1]) vertex[0] = totalAxis[0][1];
                    if(vertex[1] < totalAxis[1][0]) vertex[1] = totalAxis[1][0];
                    if(vertex[1] > totalAxis[1][1]) vertex[1] = totalAxis[1][1];
                    polygon.addPoint(new StsPoint(vertex[0], vertex[1], 0.0f), false);
                }
                polygon.close();

                // Add polygon to crossplot
                crossplot.clearCrossplotTextureDisplays();
                crossplot.polygonsChanged();
                crossplot.processPolygons();
                line = polygonFile.readLine();
                if( line == null) break;
            }
            return true;
        }
        catch(Exception e)
        {
            new StsMessage(wizard.frame, StsMessage.ERROR, "Failed to load polygons.\n" + e.getMessage());
            return false;
        }
    }
}
