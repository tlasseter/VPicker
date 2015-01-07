package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.StsSensor;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSensorRelativePanel extends StsJPanel {
    private StsSensorLoadWizard wizard;
    private StsSensorRelative wizardStep;

    StsListFieldBean fileListBean = new StsListFieldBean();
	JScrollPane fileScrollPane = new JScrollPane();
    private Object[] selectedFiles = null;

    StsGroupBox coorBox = new StsGroupBox("File Coordinates are Relative to...");
    StsDoubleFieldBean xField = new StsDoubleFieldBean();
    StsDoubleFieldBean yField = new StsDoubleFieldBean();
    StsDoubleFieldBean zField = new StsDoubleFieldBean();
    private JTextPane informationPane = new JTextPane();

    private StsSensorFile currentFile = null;
    FlowLayout flowLayout1 = new FlowLayout(FlowLayout.CENTER);


    public StsSensorRelativePanel(StsWizard wizard, StsWizardStep wizardStep) {
        this.wizard = (StsSensorLoadWizard) wizard;
        this.wizardStep = (StsSensorRelative) wizardStep;

        try {
            jbInit();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        StsProject project = wizard.getModel().getProject();
        double[] xy = project.getAbsoluteXYCoordinates(project.getXMin(), project.getYMin());

        StsSensorFile[] sensorFiles = wizard.getRelativeXYZSensorFiles();
        currentFile = sensorFiles[0];
        selectedFiles = new Object[] {currentFile};

        // If zero, set to project mins so box does not expand to accept
        for(int i=0; i<sensorFiles.length; i++)
        {
            if(sensorFiles[i].getRelativeX() == 0.0f)
                sensorFiles[i].setRelativeX(xy[0]);
            if(sensorFiles[i].getRelativeY() == 0.0f)
                sensorFiles[i].setRelativeY(xy[1]);
        }
        fileListBean.initialize(this,"files","Selected Files:",sensorFiles);
        fileListBean.setSelectedIndex(0);

        xField.initialize(this, "relativeXLoc", true, "X:");
        yField.initialize(this, "relativeYLoc", true, "Y:");
        zField.initialize(this, "relativeZLoc", true, "Z:");

        informationPane.setLayout(flowLayout1);
        informationPane.setText("Some of the selected files appear to have relative X, Y, Z values." +
                "While these can be loaded, if other project data is in absolute coordinates\n" +
                "the resulting project bounding box may be expanded and corrupt the project.\n\n" +
                "To leave the coordinates as in the file, press Next>> otherwise enter coordinates.");
        wizard.rebuild();
    }

    void jbInit() throws Exception
    {
        gbc.fill = this.gbc.HORIZONTAL;
        addEndRow(informationPane);

        fileScrollPane.getViewport().add(fileListBean, null);
        gbc.fill = this.gbc.BOTH;
        gbc.weighty = 1.0;
        addEndRow(fileScrollPane);

        coorBox.addToRow(xField);
        coorBox.addToRow(yField);
        coorBox.addEndRow(zField);

        gbc.fill = this.gbc.HORIZONTAL;
        addEndRow(coorBox);
    }

    public void setFiles(Object file)
    {
        selectedFiles = fileListBean.getSelectedObjects();
        currentFile = (StsSensorFile)selectedFiles[0];

        xField.setValue(currentFile.getRelativeX());
        yField.setValue(currentFile.getRelativeY());
        zField.setValue(currentFile.getRelativeZ());
    }
    public Object getFiles() { return currentFile; }
    public void setRelativeXLoc(double x)
    {
        if(selectedFiles == null) return;
        for(int i=0; i<selectedFiles.length; i++)
            ((StsSensorFile)selectedFiles[i]).setRelativeX(x);
    }
    public void setRelativeYLoc(double y)
    {
        if(selectedFiles == null) return;
        for(int i=0; i<selectedFiles.length; i++)
            ((StsSensorFile)selectedFiles[i]).setRelativeY(y);
    }
    public void setRelativeZLoc(double z)
    {
        if(selectedFiles == null) return;        
        for(int i=0; i<selectedFiles.length; i++)
            ((StsSensorFile)selectedFiles[i]).setRelativeZ(z);
    }
    public double getRelativeXLoc() { return currentFile.getRelativeX(); }
    public double getRelativeYLoc() { return currentFile.getRelativeY(); }
    public double getRelativeZLoc() { return currentFile.getRelativeZ(); }
}
