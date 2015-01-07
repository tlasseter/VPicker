package com.Sts.Actions.Wizards.Surfaces;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Import.*;
import com.Sts.IO.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

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

public class StsDefineRowsPanel extends StsJPanel
{
    StsJPanel parametersPanel = new StsJPanel();
    StsGroupBox fileViewPanel = new StsGroupBox("File View");
    StsComboBoxFieldBean fileComboBean;
    JTextArea fileTextArea;

    int height = 500;
    int width = 500;

    private StsAbstractFile currentFile = null;

    private StsSurfaceWizard wizard;

    int numCols = 0;
    int numSkipped = 0;

    public StsDefineRowsPanel(StsWizard wizard, StsDefineRows wizardStep)
    {
        this.wizard = (StsSurfaceWizard)wizard;

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
        StsAbstractFile[] viewableFiles = wizard.getViewableFiles();
        currentFile = wizard.getCurrentFile();
        if(!StsMath.arrayContains(viewableFiles, currentFile))
            currentFile = viewableFiles[0];
        fileComboBean.setListItems(viewableFiles, currentFile);
        updateFileView();
        repaint();
    }

    void jbInit() throws Exception
    {
        parametersPanel.gbc.fill = GridBagConstraints.HORIZONTAL;       
        fileComboBean = new StsComboBoxFieldBean(this, "file", "File:");
        parametersPanel.gbc.weightx = 1.0;
        parametersPanel.addToRow(fileComboBean);
        parametersPanel.gbc.weightx = 0.0;
        // parametersPanel.gbc.anchor = GridBagConstraints.EAST;
        StsIntFieldBean numHeaderRowsBean = new StsIntFieldBean(this, "numSkipped", 0, 1000, "Number of Header Rows", true);
        parametersPanel.addEndRow(numHeaderRowsBean);

        gbc.fill = GridBagConstraints.HORIZONTAL;
        add(parametersPanel);

        StsTextAreaScrollPane scrollPane = new StsTextAreaScrollPane();
        fileTextArea = scrollPane.textArea;
        fileViewPanel.gbc.fill = GridBagConstraints.BOTH;
        fileViewPanel.add(scrollPane);

        gbc.fill = GridBagConstraints.BOTH;
        add(fileViewPanel);
    }

    public void setNumSkipped(int n)
    {
        numSkipped = n;
        StsImportSeismicSurfaces.setSkippedLines(numSkipped);
        updateFileView();
    }

    public int getNumSkipped() { return numSkipped; }

    public void setFile(Object file)
    {
        if(!(file instanceof StsAbstractFile)) return;
        currentFile = (StsAbstractFile)file;
        updateFileView();
    }

    public Object getFile() { return currentFile; }

    public void skipHeader(StsAsciiFile asciiFile)
    {
        try
        {
            for (int i = 0; i < numSkipped; i++)
                asciiFile.readLine();
        }
        catch(Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING,  e.getMessage());
            return;
        }
    }

    public int getNumSkippedRows() { return numSkipped; }

    public boolean updateFileView()
    {
        if(currentFile == null)
            return false;

        if(!StsKeywordIO.isFileAscii(currentFile.filename))
        {
            new StsMessage(wizard.model.win3d, StsMessage.INFO,  "This is a binary file and cannot be isVisible or edited.");
            fileTextArea.setText("");
            return false;
        }
        StsAsciiFile asciiFile = new StsAsciiFile(currentFile);
        if(!asciiFile.openReadWithErrorMessage()) return false;

        try
        {
            fileTextArea.setText("");
            skipHeader(asciiFile);
            String line = asciiFile.readLine();
            int nLinesRead = 1;
            while((line != null) && (nLinesRead < 20))
            {
                appendLine(line);
                line = asciiFile.readLine();
                nLinesRead++;
            }
        }
        catch(Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "File header read error for " +
                    currentFile.filename + ": " + e.getMessage());
            return false;
        }

        return true;
    }
    public void appendLine(String line)
    {
        fileTextArea.append(line + '\n');
    }
}
