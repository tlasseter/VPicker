package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.*;
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

public class StsSensorStaticPanel extends StsJPanel {
    private StsSensorLoadWizard wizard;
    private StsSensorStatic wizardStep;

    StsListFieldBean fileListBean = new StsListFieldBean();
	JScrollPane fileScrollPane = new JScrollPane();
    private Object[] selectedFiles = null;
    
    StsGroupBox coorBox = new StsGroupBox("Specify Location of Sensor...");

    StsDoubleFieldBean xField = new StsDoubleFieldBean();
    StsDoubleFieldBean yField = new StsDoubleFieldBean();
    StsDoubleFieldBean zField = new StsDoubleFieldBean();
    private JTextPane informationPane = new JTextPane();

    private StsSensorFile currentFile = null;
    FlowLayout flowLayout1 = new FlowLayout(FlowLayout.CENTER);

    public StsSensorStaticPanel(StsWizard wizard, StsWizardStep wizardStep) {
        this.wizard = (StsSensorLoadWizard) wizard;
        this.wizardStep = (StsSensorStatic) wizardStep;

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

        StsSensorFile[] sensorFiles = wizard.getStaticSensorFiles();
        currentFile = sensorFiles[0];
        selectedFiles = new Object[] {currentFile};

        // If zero, set to project mins so box does not expand to accept
        for(int i=0; i<sensorFiles.length; i++)
        {
            if(sensorFiles[i].getStaticX() == 0.0f)
                sensorFiles[i].setStaticX(xy[0]);
            if(sensorFiles[i].getStaticY() == 0.0f)
                sensorFiles[i].setStaticY(xy[1]);
        }
        fileListBean.initialize(this,"files","Selected Files:",sensorFiles);
        fileListBean.setSelectedIndex(0);

        xField.initialize(this, "staticXLoc", true, "X:");
        yField.initialize(this, "staticYLoc", true, "Y:");
        zField.initialize(this, "staticZLoc", true, "Z:");

        informationPane.setLayout(flowLayout1);
        informationPane.setText("Some of the selected files do not contain X, Y, Z values." +
                "A static position will need to be supplied.\n\n" +
                "We are assuming all selected files without positions are located at the same place." +
                "If not, return to the file selection step and only select files with the same location.");
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

        xField.setValue(currentFile.getStaticX());
        yField.setValue(currentFile.getStaticY());
        zField.setValue(currentFile.getStaticZ());
    }
    
    public Object getFiles() { return currentFile; }
    public void setStaticXLoc(double x)
    {
        if(selectedFiles == null) return;
        for(int i=0; i<selectedFiles.length; i++)
            ((StsSensorFile)selectedFiles[i]).setStaticX(x);
    }
    public void setStaticYLoc(double y)
    {
        if(selectedFiles == null) return;
        for(int i=0; i<selectedFiles.length; i++)
            ((StsSensorFile)selectedFiles[i]).setStaticY(y);
    }
    public void setStaticZLoc(double z)
    {
        if(selectedFiles == null) return;        
        for(int i=0; i<selectedFiles.length; i++)
            ((StsSensorFile)selectedFiles[i]).setStaticZ(z);
    }
    public double getStaticXLoc() { return currentFile.getStaticX(); }
    public double getStaticYLoc() { return currentFile.getStaticY(); }
    public double getStaticZLoc() { return currentFile.getStaticZ(); }

}
