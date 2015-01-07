package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.DBTypes.StsColor;
import com.Sts.IO.StsAbstractFile;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;
import java.util.Date;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSensorDefineColorsPanel extends StsJPanel {
    private StsSensorLoadWizard wizard;
    private StsSensorDefineColors wizardStep;

    private StsComboBoxFieldBean fileComboBean;
    StsSensorFile currentFile = null;

    //StsColorListFieldBean[] colorListBeans = null;
	JScrollPane scrollPane = new JScrollPane();
    StsGroupBox stageBox = new StsGroupBox("Define Stages");

    private StsColor[] colors = null;
    private long[] startTimes = null;
    private long[] endTimes = null;
    StsDateFieldBean[] startBeans = null;
    StsDateFieldBean[] endBeans = null;
    StsIntFieldBean numStagesBean = new StsIntFieldBean();
    int numStages = 1;

    long oneDay = (long)(1000*60*60*24);

    public StsSensorDefineColorsPanel(StsWizard wizard, StsWizardStep wizardStep) {
        this.wizard = (StsSensorLoadWizard) wizard;
        this.wizardStep = (StsSensorDefineColors) wizardStep;

        try {
            jbInit();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        StsSensorFile[] files = wizard.getSensorFiles();
        currentFile = files[0];
        fileComboBean.setListItems(files, currentFile);

        configureBeans();
    }

    private void configureBeans()
    {
        stageBox.removeAll();

        numStages = currentFile.numStages(wizard.getModel());
        numStagesBean.setValue(numStages);

        startBeans = new StsDateFieldBean[numStages];
        endBeans = new StsDateFieldBean[numStages];
        //colorListBeans = new StsColorListFieldBean[numStages];

        for(int i=0; i<numStages; i++)
        {
            startBeans[i] = new StsDateFieldBean(this, "startTimeString", StsDateFieldBean.convertToString(currentFile.getStartTimeForStage(wizard.getModel(), i)), true, "Start Time:");
            endBeans[i] = new StsDateFieldBean(this, "endTimeString", StsDateFieldBean.convertToString(currentFile.getEndTimeForStage(wizard.getModel(), i)), true, "End Time:");
            //colorListBeans[i] = new StsColorListFieldBean(this, "color", "Color: ", StsColor.colors32);
            stageBox.addToRow(startBeans[i]);
            stageBox.addToRow(endBeans[i]);
            //stageBox.addEndRow(colorListBeans[i]);
        }
    }

    void jbInit() throws Exception
    {
        gbc.fill = gbc.HORIZONTAL;
        fileComboBean = new StsComboBoxFieldBean(this, "file", "File:");
        numStagesBean.initialize(this,"numStages", 1, 100, "Number of Stages:", true);
        addEndRow(fileComboBean);
        addEndRow(numStagesBean);

        scrollPane.getViewport().add(stageBox, null);
        gbc.fill = this.gbc.BOTH;
        gbc.weighty = 1.0;
        addEndRow(scrollPane);
    }

    public StsColor getColor()
    {
        return null;
    }
    public void setColor(StsColor color)
    {
    }

    public String getStartTimeString()
    {
        return null;
    }
    public void setStartTimeString(String time)
    {
    }


    public void setFile(Object file)
    {
        if(!(file instanceof StsSensorFile)) return;
        currentFile = (StsSensorFile)file;
        configureBeans();
    }

    public StsSensorFile getFile() { return currentFile; }

    public String getEndTimeString()
    {
        return null;
    }
    public void setEndTimeString(String time)
    {
        // ToDo: Find the line in the file that corresponds to the specified time.

    }

    public int getNumStages()
    {
        return numStages;
    }
    public void setNumStages(int num)
    {
        numStages = num;
        // Reconfigure beans and rows
        currentFile.setNumStages(numStages);
        configureBeans();
    }
}