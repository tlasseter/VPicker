package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.DBTypes.StsColor;
import com.Sts.IO.StsAsciiFile;
import com.Sts.UI.Beans.*;
import com.Sts.UI.StsMessage;
import com.Sts.Utilities.StsMath;

import javax.swing.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSensorDefineStagesPanel extends StsJPanel
{
    private StsSensorLoadWizard wizard;
    private StsSensorDefineStages wizardStep;

    //private StsComboBoxFieldBean fileComboBean;
    JScrollPane fileScrollPane = new JScrollPane();
    StsListFieldBean fileListBean;
    StsSensorFile currentFile = null;
    int[] stageOffsets = new int[1];
	JScrollPane scrollPane = new JScrollPane();
    JScrollPane filePane = new JScrollPane();
    JTextArea fileTextArea = new JTextArea();

    StsGroupBox stageBox = new StsGroupBox();

    private int[] startRows = null;
    StsIntFieldBean[] startBeans = null;
    StsIntFieldBean numStagesBean = new StsIntFieldBean();
    int numStages = 1;

    public StsSensorDefineStagesPanel(StsWizard wizard, StsWizardStep wizardStep) {
        this.wizard = (StsSensorLoadWizard) wizard;
        this.wizardStep = (StsSensorDefineStages) wizardStep;

        try
        {
            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        StsSensorFile[] files = wizard.getSensorFiles();
        for(int i=0; i<files.length; i++)
        {
            files[i].numStages(wizard.getModel());
            files[i].setOverrideStages(true);
        }
        currentFile = files[0];
        fileListBean.initialize(this,"file",null,files);
        fileListBean.setSingleSelect();
        fileListBean.setSelectedValue(currentFile);
        configureBeans();
    }

    private void configureBeans()
    {
        stageBox.removeAll();

        numStages = currentFile.numStages(wizard.getModel());
        stageOffsets = currentFile.getStageOffsets();

        numStagesBean.setValue(numStages);
        stageBox.addEndRow(numStagesBean);

        startBeans = new StsIntFieldBean[numStages];
        for(int i=0; i<numStages; i++)
        {
            startBeans[i] = new StsIntFieldBean(this, "startRow"+i, 1, 1000, "Start Row Stage #:" + (i+1), true);
            stageBox.addEndRow(startBeans[i]);
        }
        fileTextFill();
        wizard.rebuild();
    }

    private boolean fileTextFill()
    {
        if(currentFile == null)
            return false;

        StsAsciiFile asciiFile = new StsAsciiFile(currentFile.file);
        if(!asciiFile.openReadWithErrorMessage())
            return false;

        try
        {
            fileTextArea.setText("");
            String line = asciiFile.readLine();
            int nLinesRead = 1;
            while(line != null)
            {
                fileTextArea.append(nLinesRead + " --> " + line + "\n");
                line = asciiFile.readLine();
                nLinesRead++;
            }
        }
        catch(Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "File header read error for " +
                    currentFile.file.filename + ": " + e.getMessage());
            return false;
        }
        return true;
    }

    void jbInit() throws Exception
    {
        fileListBean = new StsListFieldBean();
        JScrollPane fileScrollPane = new JScrollPane();
        fileScrollPane.getViewport().add(fileListBean, null);

        StsJPanel parametersPanel = new StsGroupBox("Select File to View");
        parametersPanel.gbc.fill = gbc.BOTH;
        parametersPanel.gbc.weighty = 1.0;
        parametersPanel.addEndRow(fileScrollPane);
        parametersPanel.setPreferredSize(300, 150);

        gbc.anchor = gbc.NORTH;
        gbc.weighty = 0.25;
        gbc.fill = gbc.HORIZONTAL;
        addEndRow(parametersPanel);

        StsJPanel defineStagesPanel = new StsGroupBox("Define Stages");
        filePane.getViewport().add(fileTextArea, null);
        scrollPane.getViewport().add(stageBox, null);

        numStagesBean.initialize(this,"numStages", 1, 20, "Number of Stages:", true);
        stageBox.addEndRow(numStagesBean);
        
        defineStagesPanel.gbc.fill = this.gbc.BOTH;
        defineStagesPanel.gbc.weighty = 1.0;
        defineStagesPanel.addToRow(scrollPane);
        defineStagesPanel.addEndRow(filePane);

        defineStagesPanel.setPreferredSize(300, 300);        
        gbc.fill = gbc.BOTH;
        gbc.weighty = 0.75;
        addEndRow(defineStagesPanel);
    }

    public int getStartRow0() { return stageOffsets[0] + 1; }
    public void setStartRow0(int row) { stageOffsets[0] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow1() { return stageOffsets[1] + 1; }
    public void setStartRow1(int row) { stageOffsets[1] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow2() { return stageOffsets[2] + 1; }
    public void setStartRow2(int row) { stageOffsets[2] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow3() { return stageOffsets[3] + 1; }
    public void setStartRow3(int row) { stageOffsets[3] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow4() { return stageOffsets[4] + 1; }
    public void setStartRow4(int row) { stageOffsets[4] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow5() { return stageOffsets[5] + 1; }
    public void setStartRow5(int row) { stageOffsets[5] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow6() { return stageOffsets[6] + 1; }
    public void setStartRow6(int row) { stageOffsets[6] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow7() { return stageOffsets[7] + 1; }
    public void setStartRow7(int row) { stageOffsets[7] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow8() { return stageOffsets[8] + 1; }
    public void setStartRow8(int row) { stageOffsets[8] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow9() { return stageOffsets[9] + 1; }
    public void setStartRow9(int row) { stageOffsets[9] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow10() { return stageOffsets[10] + 1; }
    public void setStartRow10(int row) { stageOffsets[10] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow11() { return stageOffsets[11] + 1; }
    public void setStartRow11(int row) { stageOffsets[11] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow12() { return stageOffsets[12] + 1; }
    public void setStartRow12(int row) { stageOffsets[12] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow13() { return stageOffsets[13] + 1; }
    public void setStartRow13(int row) { stageOffsets[13] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow14() { return stageOffsets[14] + 1; }
    public void setStartRow14(int row) { stageOffsets[14] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow15() { return stageOffsets[15] + 1; }
    public void setStartRow15(int row) { stageOffsets[15] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow16() { return stageOffsets[16] + 1; }
    public void setStartRow16(int row) { stageOffsets[16] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow17() { return stageOffsets[17] + 1; }
    public void setStartRow17(int row) { stageOffsets[17] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow18() { return stageOffsets[18] + 1; }
    public void setStartRow18(int row) { stageOffsets[18] = row - 1; currentFile.setStageOffsets(stageOffsets); }
    public int getStartRow19() { return stageOffsets[19] + 1; }
    public void setStartRow19(int row) { stageOffsets[19] = row - 1; currentFile.setStageOffsets(stageOffsets); }

    public void setFile(StsSensorFile file)
    {
        if(!(file instanceof StsSensorFile)) return;
        currentFile = (StsSensorFile)file;
        configureBeans();
    }

    public StsSensorFile getFile() { return currentFile; }

    public int getNumStages()
    {
        return numStages;
    }

    public void setNumStages(int num)
    {
        if(num == numStages) return;
        int diff = num - numStages;
        int eof = stageOffsets[numStages];
        if(diff > 0)
        {
            stageOffsets = (int[]) StsMath.arrayGrow(stageOffsets, diff);
            for(int i=numStages; i<=num; i++)
                stageOffsets[i] = 1;
        }
        else
        {
            stageOffsets = (int[]) StsMath.arrayDeleteElementRange(stageOffsets, num, numStages);
        }
        stageOffsets[num] = eof;
        numStages = num;

        // Reconfigure beans and rows
        currentFile.setNumStages(numStages);
        currentFile.setStageOffsets(stageOffsets);
        configureBeans();
    }
}