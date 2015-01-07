package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.Beans.StsDateFieldBean;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Beans.StsListFieldBean;

import javax.swing.*;
import java.awt.*;
import java.text.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSensorTimePanel extends StsJPanel {
    private StsSensorLoadWizard wizard;
    private StsSensorTime wizardStep;

    StsListFieldBean fileListBean = new StsListFieldBean();
	JScrollPane fileScrollPane = new JScrollPane();

    private StsSensorFile currentFile = null;
    private Object[] selectedFiles = null;
    StsGroupBox coorBox = new StsGroupBox("Times in file are relative to...");
    StsDateFieldBean dateBean = new StsDateFieldBean();

    private JTextPane informationPane = new JTextPane();
    FlowLayout flowLayout1 = new FlowLayout(FlowLayout.CENTER);
    
    long oneDay = (long)(1000*60*60*24);
                                 
    public StsSensorTimePanel(StsWizard wizard, StsWizardStep wizardStep) {
        this.wizard = (StsSensorLoadWizard) wizard;
        this.wizardStep = (StsSensorTime) wizardStep;

        try {
            jbInit();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        StsSensorFile[] viewableFiles = wizard.getRelativeTimeSensorFiles();
        currentFile = viewableFiles[0];
        fileListBean.initialize(this,"files","Selected Files:",viewableFiles);
        fileListBean.setSelectedIndex(0);
        selectedFiles = new Object[] {currentFile};

        informationPane.setLayout(flowLayout1);
        informationPane.setText("Some of the selected files only contain relative times, therefore date (dd-mm-yy) is required.\n\n" +
                "We are assuming all selected files without a date are from the same day." +
                " If not, return to the file selection step and only select sensors collected on the same day.");
        
        Date date = new Date(oneDay);
        dateBean.setFormat(wizard.model.getProject().getTimeDateFormat());
        String dateString = dateBean.convertToString(date.getTime());

        dateBean.initialize(this,"startTimeString",dateString,true,"Absolute Start Time:");
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

        coorBox.addEndRow(dateBean);
        gbc.fill = this.gbc.HORIZONTAL;
        addEndRow(coorBox);
    }

    public String getStartTimeString()
    {
        return StsDateFieldBean.convertToString(currentFile.getStartTime());
    }
    public void setStartTimeString(String time)
    {
        if(selectedFiles == null) return;
        for(int i=0; i<selectedFiles.length; i++)
        {
            ((StsSensorFile)selectedFiles[i]).setStartTime(dateBean.convertToLong(time));
        }
    }

    public void setFiles(Object file)
    {
        selectedFiles = fileListBean.getSelectedObjects();
        currentFile = (StsSensorFile)selectedFiles[0];
        dateBean.setValue(dateBean.convertToString(currentFile.getStartTime()));
    }

    public Object getFiles() { return currentFile; }
}
