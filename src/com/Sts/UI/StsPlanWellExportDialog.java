package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsPlanWellExportDialog extends StsWellExportDialog
{
    private boolean exportMarkers = false;

    StsGroupBox markersBox;
	StsBooleanFieldBean exportMarkersBean;
    ArrayList markerObjectsList;
    StsTablePanelNew markerTable;
    StsAsciiFile asciiFile = null;

    public StsPlanWellExportDialog(StsModel model, Frame frame, String title, boolean modal, StsWell well)
    {
        super(model, frame, title, modal, well, "");
    }


    static public boolean exportWell(StsModel model, Frame frame, String title, boolean modal, StsWell well, String timeOrDepth)
    {
        try
        {
            StsPlanWellExportDialog dialog = new StsPlanWellExportDialog(model, frame, title, modal, well);
            dialog.setSize(200, 600);
            dialog.pack();
            dialog.setVisible(true);
            dialog.exportWell(timeOrDepth);
            dialog.exportHorizonMarkers();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsWellExportDialog.class, "constructor", e);
            return false;
        }
    }

    protected void constructBeans()
    {
        super.constructBeans();
        markersBox = new StsGroupBox("Export markers");
	    exportMarkersBean = new StsBooleanFieldBean(this, "exportMarkers", "Export markers");
    }

    protected void constructPanel()
    {
        constructMarkerObjects();       
        this.getContentPane().add(panel);
		this.setTitle("Well Export Parameters");
		panel.add(selectBox);
        panel.add(markersBox);
        panel.add(rangeBox);
		panel.add(buttonPanel);

		selectBox.add(nameBean);
		selectBox.addToRow(logsDataBean);
		selectBox.addToRow(deviationDataBean);
		selectBox.addEndRow(seisAttsDataBean);

        markersBox.add(exportMarkersBean);
        markersBox.add(markerTable);

        rangeBox.addToRow(minMDepthBean);
		rangeBox.addEndRow(maxMDepthBean);

	    buttonPanel.addToRow(processButton);
		buttonPanel.addEndRow(cancelButton);

	    minMDepthBean.setValue(0.0f);
        minMDepth = 0.0f;
		maxMDepthBean.setValue(well.getLastPoint().getM());
        maxMDepth = well.getLastPoint().getM();
    }

    /** exports a well-ref.txt.wellname markers file.  Returns true if no markers exist or executes successfully. */
    private boolean exportHorizonMarkers()
    {
        if(!exportMarkers) return true;
        int nMarkers = markerObjectsList.size();
        if(nMarkers == 0) return true;
        String filename = "well-ref.txt." + exportName;
        String pathname = model.getProject().getRootDirString() + filename;

        try
        {
           StsFile stsFile = StsFile.constructor(pathname);
           asciiFile = new StsAsciiFile(stsFile);

            if(!asciiFile.openWrite())
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "Failed to open file for writing: " + pathname);
                return false;
            }
            asciiFile.writeLine("WELLNAME");
            asciiFile.writeLine(exportName);
            asciiFile.writeLine("CURVE");
            asciiFile.writeLine("MDEPTH");
            asciiFile.writeLine("VALUE");
            writeMarkerObjects();
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "exportMarkers", e);
            return false;
        }
        finally
        {
            if(asciiFile != null) asciiFile.close();
        }
    }


    private void constructMarkerObjects()
    {
        StsMainObject[] markers = model.getObjectListOfType(StsMarker.class, StsMarker.SURFACE);
        int nMarkers = markers.length;
        markerObjectsList = new ArrayList();
        for(int n = 0; n < nMarkers; n++)
        {
            StsMarker marker = (StsMarker)markers[n];
            StsModelSurface modelSurface = marker.getModelSurface();
            float intersection = StsWellMarker.computeWellSurfaceIntersection(modelSurface, well);
            markerObjectsList.add(new MarkerTableObject(marker.getName(), intersection));
        }
        String[] colNames = new String[] { "name", "depth"};
        String[] colTitles = new String[] { "Marker name", "MDepth"};
        markerTable = new StsTablePanelNew(markerObjectsList, colNames, colTitles);
        markerTable.addDeleteButtons();
    }

    private void writeMarkerObjects()
    {
        try
        {
            int nMarkers = markerObjectsList.size();
            for(int n = 0; n < nMarkers; n++)
            {
                MarkerTableObject markerObject = (MarkerTableObject)markerObjectsList.get(n);
//                if(markerObject.depth != StsParameters.nullValue)
                    asciiFile.writeLine(markerObject.name + " " + Float.toString(markerObject.depth));
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "writeMarkerObjects", e);
        }
    }

    public boolean getExportMarkers() { return exportMarkers; }
    public void setExportMarkers(boolean exportMarkers) { this.exportMarkers = exportMarkers; }

    public class MarkerTableObject
    {
        public String name;
        public float depth;

        public MarkerTableObject(String name, float depth)
        {
            this.name = name;
            this.depth = depth;
        }
    }
}