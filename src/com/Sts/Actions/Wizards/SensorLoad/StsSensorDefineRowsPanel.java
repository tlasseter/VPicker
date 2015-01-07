package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Import.StsKeywordIO;
import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.IO.StsAbstractFile;
import com.Sts.IO.StsAsciiFile;
import com.Sts.UI.Beans.*;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsTextAreaScrollPane;
import com.Sts.Utilities.StsMath;

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

public class StsSensorDefineRowsPanel extends StsJPanel
{

    JTextArea fileTextArea;
    StsIntFieldBean numHeaderRowsBean;
    StsListFieldBean fileListBean;
    private Object[] selectedFiles = null;

    private StsSensorFile currentFile = null;

    private StsSensorLoadWizard wizard;

    int numCols = 0;

    public StsSensorDefineRowsPanel(StsWizard wizard, StsSensorDefineRows wizardStep)
    {
        this.wizard = (StsSensorLoadWizard)wizard;

        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        StsSensorFile[] viewableFiles = wizard.getSensorFiles();
        currentFile = viewableFiles[0];
        selectedFiles = new Object[] {currentFile};

        numHeaderRowsBean.getValueFromPanelObject();
        fileListBean.initialize(this,"files",null,viewableFiles);
        fileListBean.setSelectedIndex(0);
        updateFileView();
        repaint();

        fileListBean.setSelectedAll();
    }

    void jbInit() throws Exception
    {
        fileListBean = new StsListFieldBean();
        JScrollPane fileScrollPane = new JScrollPane();
        fileScrollPane.getViewport().add(fileListBean, null);

        StsJPanel parametersPanel = new StsGroupBox("Select File to View");
        parametersPanel.setPreferredSize(300, 150);
        parametersPanel.gbc.fill = gbc.BOTH;
        parametersPanel.gbc.weighty = 1.0;
        parametersPanel.addEndRow(fileScrollPane);

        numHeaderRowsBean = new StsIntFieldBean(this, "numSkipped", 0, 1000, "Number of Header Rows", true);
        parametersPanel.gbc.fill = gbc.HORIZONTAL;
        parametersPanel.gbc.weighty = 0.0;
        parametersPanel.addEndRow(numHeaderRowsBean);

        gbc.fill = GridBagConstraints.BOTH;
        gbc.weighty = 0.8;
        add(parametersPanel);

        StsTextAreaScrollPane scrollPane = new StsTextAreaScrollPane();
        fileTextArea = scrollPane.textArea;
        fileTextArea.setPreferredSize(new Dimension(300, 150));
        StsGroupBox fileViewPanel = new StsGroupBox("File View");
        fileViewPanel.gbc.fill = GridBagConstraints.BOTH;
        fileViewPanel.add(scrollPane);
        gbc.weighty = 0.2;
        add(fileViewPanel);
    }

    public void setNumSkipped(int n)
    {
        if(selectedFiles == null) return;        
        for(int i=0; i<selectedFiles.length; i++)
           ((StsSensorFile)selectedFiles[i]).setNumHeaderRows(n);
        updateFileView();
    }

    public int getNumSkipped() { return currentFile.getNumHeaderRows(); }

    public void setFiles(Object file)
    {
        selectedFiles = fileListBean.getSelectedObjects();
        currentFile = (StsSensorFile)selectedFiles[0];
        numHeaderRowsBean.setValue(currentFile.getNumHeaderRows());
        updateFileView();
    }

    public Object getFiles() { return currentFile; }

    public void skipHeader(StsAsciiFile asciiFile)
    {
        try
        {
            for (int i = 0; i < currentFile.getNumHeaderRows(); i++)
                asciiFile.readLine();
        }
        catch(Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING,  e.getMessage());
            return;
        }
    }

    public boolean updateFileView()
    {
        if(currentFile == null)
            return false;

        StsAsciiFile asciiFile = new StsAsciiFile(currentFile.file);
        if(!asciiFile.openReadWithErrorMessage()) return false;

        try
        {
            fileTextArea.setText("");
            skipHeader(asciiFile);
            String line = asciiFile.readLine();
            int nLinesRead = 1;
            while((line != null) && (nLinesRead < 15))
            {
                appendLine(line);
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
    public void appendLine(String line)
    {
        fileTextArea.append(line + '\n');
    }
}