package com.Sts.Actions.Wizards.Surfaces;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Import.*;
import com.Sts.IO.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineColumnsPanel extends StsJPanel implements ActionListener, ChangeListener
{
    private StsTablePanel fileTable = new StsTablePanel();
	private JScrollPane tableScrollPane1 = new JScrollPane();
    private JLabel numColLabel = new JLabel("# Cols:");
    // private JComboBox fileCombo = new JComboBox();
    private StsComboBoxFieldBean fileComboBean;

    private StsGroupBox colDefPanel = new StsGroupBox("Column Assign");
    private StsAbstractFile currentFile = null;
    private StsSurfaceWizard wizard;
    private StsDefineColumns wizardStep;

    int selectedCol, xCol, yCol, tCol, ilineCol, xlineCol;
    int numCols = 0;
    int numSkipped = 0;

    byte X = StsImportSeismicSurfaces.X;
    byte Y = StsImportSeismicSurfaces.Y;
    byte T = StsImportSeismicSurfaces.T;
    byte ILINE = StsImportSeismicSurfaces.ILINE;
    byte XLINE = StsImportSeismicSurfaces.XLINE;


    //JPanel colDefPanel = new JPanel();
    JLabel jLabel2 = new JLabel();
    JSpinner numColSpin = new JSpinner();
    private SpinnerModel numColSpinModel = null;
    JSpinner xSpin = new JSpinner();
    private SpinnerModel xSpinModel = null;
    JSpinner ySpin = new JSpinner();
    private SpinnerModel ySpinModel = null;
    JSpinner tSpin = new JSpinner();
    private SpinnerModel tSpinModel = null;
    JSpinner inlineSpin = new JSpinner();
    private SpinnerModel inlineSpinModel = null;
    JSpinner xlineSpin = new JSpinner();
    private SpinnerModel xlineSpinModel = null;
    JLabel jLabel3 = new JLabel();
    JLabel jLabel7 = new JLabel();
    JLabel jLabel8 = new JLabel();
    JLabel jLabel9 = new JLabel();
    JLabel jLabel10 = new JLabel();
    JLabel jLabel11 = new JLabel();
    JLabel jLabel12 = new JLabel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JCheckBox concatLinesChk = new JCheckBox("Concatenate Short Lines");
    boolean concatLines = true;
    
    StsGroupBox fileBox = new StsGroupBox();
    //JPanel jPanel1 = new JPanel();
    JPanel jPanel2 = new JPanel();
    JButton okBtn = new JButton();
    JButton cancelBtn = new JButton();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    GridBagLayout gridBagLayout4 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsDefineColumnsPanel(StsWizard wizard, StsDefineColumns wizardStep)
    {
        this.wizard = (StsSurfaceWizard)wizard;
        this.wizardStep = wizardStep;

        try
        {
            constructPanel();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    void constructPanel() throws Exception
    {
        this.setLayout(gridBagLayout4);

        numColSpinModel = new SpinnerNumberModel(StsImportSeismicSurfaces.numColsPerRow, 1, 200, 1);
        xSpinModel = new SpinnerNumberModel(StsImportSeismicSurfaces.colOrder[0]+1, 1, 200, 1);
        ySpinModel = new SpinnerNumberModel(StsImportSeismicSurfaces.colOrder[1]+1, 1, 200, 1);
        tSpinModel = new SpinnerNumberModel(StsImportSeismicSurfaces.colOrder[2]+1, 1, 200, 1);
        inlineSpinModel = new SpinnerNumberModel(StsImportSeismicSurfaces.colOrder[3]+1, 1, 200, 1);
        xlineSpinModel = new SpinnerNumberModel(StsImportSeismicSurfaces.colOrder[4]+1, 1, 200, 1);

        fileComboBean = new StsComboBoxFieldBean(this, "file", "File:");

        fileTable.setTitle("Surface File:");
		tableScrollPane1.setAutoscrolls(true);
		tableScrollPane1.getViewport().add(fileTable, null);
		tableScrollPane1.getViewport().setPreferredSize(new Dimension(550, 100));
		
        //colDefPanel.setBorder(BorderFactory.createEtchedBorder());
        //colDefPanel.setLayout(gridBagLayout1);
        jLabel3.setFont(new java.awt.Font("Serif", 1, 11));
        jLabel3.setHorizontalAlignment(SwingConstants.CENTER);
        jLabel3.setText("Field");
        jLabel7.setFont(new java.awt.Font("Serif", 1, 11));
        jLabel7.setHorizontalAlignment(SwingConstants.CENTER);
        jLabel7.setText("Column");
        jLabel8.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel8.setText("X:");
        jLabel9.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel9.setText("Y:");
        jLabel10.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel10.setText("T:");
        jLabel11.setText("Inline:");
        jLabel12.setText("Crossline:");
        numColSpin.addChangeListener(this);
        numColSpin.setModel(numColSpinModel);
        xSpin.addChangeListener(this);
        xSpin.setModel(xSpinModel);
        ySpin.addChangeListener(this);
        ySpin.setModel(ySpinModel);
        tSpin.addChangeListener(this);
        tSpin.setModel(tSpinModel);
        inlineSpin.addChangeListener(this);
        inlineSpin.setModel(inlineSpinModel);
        xlineSpin.addChangeListener(this);
        xlineSpin.setModel(xlineSpinModel);
        concatLinesChk.addActionListener(this);
        //jPanel1.setBorder(BorderFactory.createEtchedBorder());
        //jPanel1.setLayout(gridBagLayout2);
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout3);
        okBtn.setText("Ok");
        okBtn.addActionListener(this);
        cancelBtn.setText("Cancel");
        cancelBtn.addActionListener(this);

        colDefPanel.addToRow(numColLabel);
        colDefPanel.addEndRow(numColSpin);
        colDefPanel.gbc.gridwidth = 2;
        colDefPanel.addEndRow(concatLinesChk);
        colDefPanel.gbc.gridwidth = 1;
        colDefPanel.addToRow(jLabel8);
        colDefPanel.addEndRow(xSpin);
        colDefPanel.addToRow(jLabel9);
        colDefPanel.addEndRow(ySpin);
        colDefPanel.addToRow(jLabel10);
        colDefPanel.addEndRow(tSpin);
        colDefPanel.addToRow(jLabel11);
        colDefPanel.addEndRow(inlineSpin);
        colDefPanel.addToRow(jLabel12);
        colDefPanel.addEndRow(xlineSpin);

        fileBox.gbc.weighty = 0.0;
        fileBox.gbc.anchor = fileBox.gbc.WEST;
        fileBox.gbc.fill = fileBox.gbc.NONE;
        fileBox.gbc.weightx = 1.0;
        fileBox.gbc.fill = fileBox.gbc.HORIZONTAL;
        fileBox.gbc.gridwidth = 2;
        fileBox.addEndRow(fileComboBean);
        fileBox.gbc.weighty = 1.0;
        fileBox.gbc.weightx = 0.0;
        fileBox.gbc.fill = fileBox.gbc.VERTICAL;
        fileBox.gbc.gridwidth = 1;
        fileBox.addToRow(colDefPanel);
        fileBox.gbc.weightx = 1.0;
        fileBox.gbc.fill = fileBox.gbc.BOTH;
        fileBox.addEndRow(tableScrollPane1);

    this.add(fileBox,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
        ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
    }

    public void initialize()
    {
        StsImportSeismicSurfaces.gridType = StsImportSeismicSurfaces.XYZIX; // Default grid type unless header is detected that counters
        numSkipped = wizard.defineRows.panel.getNumSkippedRows();

        StsAbstractFile[] viewableFiles = wizard.getViewableFiles();
        currentFile = wizard.getCurrentFile();
        if(!StsMath.arrayContains(viewableFiles, currentFile))
            currentFile = viewableFiles[0];
        fileComboBean.setListItems(viewableFiles, currentFile);

        initNumCols();
        numColSpin.setValue(new Integer(numCols));

        xCol = StsImportSeismicSurfaces.colOrder[0];
        yCol = StsImportSeismicSurfaces.colOrder[1];
        tCol = StsImportSeismicSurfaces.colOrder[2];
        ilineCol = StsImportSeismicSurfaces.colOrder[3];
        xlineCol = StsImportSeismicSurfaces.colOrder[4];

        xSpin.setValue(new Integer(xCol+1));
        ySpin.setValue(new Integer(yCol+1));
        tSpin.setValue(new Integer(tCol+1));
        inlineSpin.setValue(new Integer(ilineCol+1));
        xlineSpin.setValue(new Integer(xlineCol+1));

        if(StsImportSeismicSurfaces.gridType == StsImportSeismicSurfaces.XYZ)
        {
            inlineSpin.setEnabled(false);
            xlineSpin.setEnabled(false);
        }
        concatLinesChk.setSelected(concatLines);
        updateFileView();
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == concatLinesChk)
        {
        	concatLines = concatLinesChk.isSelected();
        }
    }

    public void stateChanged(ChangeEvent e)
    {
        Object source = e.getSource();
        if(source == xSpin)
        {
            xCol = Integer.valueOf(xSpin.getValue().toString()).intValue() - 1;
            StsImportSeismicSurfaces.setColLocation(X, xCol);
//                wizard.setOrder(StsImportSeismicSurfaces.X, xCol);
        }
        else if(source == ySpin)
        {
            yCol = Integer.valueOf(ySpin.getValue().toString()).intValue() - 1;
            StsImportSeismicSurfaces.setColLocation(Y, yCol);
//                wizard.setOrder(StsImportSeismicSurfaces.Y, yCol);
        }
        else if(source == tSpin)
        {
            tCol = Integer.valueOf(tSpin.getValue().toString()).intValue() - 1;
            StsImportSeismicSurfaces.setColLocation(T, tCol);
//                 wizard.setOrder(StsImportSeismicSurfaces.T, tCol);
        }
        else if(source == inlineSpin)
        {
            ilineCol = Integer.valueOf(inlineSpin.getValue().toString()).intValue() - 1;
            StsImportSeismicSurfaces.setColLocation(ILINE, ilineCol);
//                 wizard.setOrder(StsImportSeismicSurfaces.ILINE, ilineCol);
        }
        else if(source == xlineSpin)
        {
            xlineCol = Integer.valueOf(xlineSpin.getValue().toString()).intValue() - 1;
            StsImportSeismicSurfaces.setColLocation(XLINE, xlineCol);
//                wizard.setOrder(StsImportSeismicSurfaces.XLINE, xlineCol);
        }
        else if(source == numColSpin)
        {
            numCols = Integer.valueOf(numColSpin.getValue().toString()).intValue();
            StsImportSeismicSurfaces.setNumCols(numCols);
        }
        if(numCols > 0)
            updateFileView();
    }

    public void setFile(Object file)
    {
        if(!(file instanceof StsAbstractFile)) return;
        currentFile = (StsAbstractFile)file;
        updateFileView();
    }

    public Object getFile() { return currentFile; }
    public boolean initNumCols()
    {
        String[] tokens = null;

        if(!StsKeywordIO.isFileAscii(currentFile.filename)) return false;
        StsAsciiFile asciiFile = new StsAsciiFile(currentFile);
        if (!asciiFile.openReadWithErrorMessage())return false;

        try
        {
            skipHeader(asciiFile, true);
            tokens = asciiFile.getTokens(new String[] {" ",","});
            numCols = tokens.length;
        }
        catch(Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "File header read error for " +
                           currentFile.filename + ": " + e.getMessage());
            return false;
        }
        return true;
    }

    public void skipHeader(StsAsciiFile asciiFile, boolean determineOrder)
    {
        String line = null;
        String timeOrDepth = wizard.selectProcessSurfaces.panel.verticalUnitsPanel.timeDepthString;
        boolean isTrueDepth = timeOrDepth.equals(StsParameters.TD_DEPTH_STRING);
        try
        {
            for (int i = 0; i < numSkipped; i++)
            {
            	// Assume last skipped line is the column names
                line = asciiFile.readLine();
                if((i == (numSkipped - 1)) && determineOrder)
                	StsImportSeismicSurfaces.determineColumnOrder(line, isTrueDepth);
            }
        }
        catch(Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING,  e.getMessage());
            return;
        }
    }

    public boolean concatLines() { return concatLines; }
    
    public boolean updateFileView()
    {
        String[] tokens;

        if(!StsKeywordIO.isFileAscii(currentFile.filename)) return false;
        StsAsciiFile asciiFile = new StsAsciiFile(currentFile);
        if(!asciiFile.openReadWithErrorMessage()) return false;

        try
        {
            skipHeader(asciiFile, false);
            tokens = asciiFile.getTokens(new String[] {" ",","});
            int nLines = 0;
            Object[] row = new Object[numCols];

            fileTable.removeAllRows();
            setColumnHeadings();

            int rowCount = 0;
            while(tokens != null)
            {
                if(++nLines > 100)
                    break;
                /*
                if (nTokens < 5)
                {
                    String inputLine = asciiFile.getLine();
                    StsMessageFiles.errorMessage("Insufficient entries for line: " + inputLine);
                    tokens = asciiFile.getTokens();
                    continue;
                }
*/
                while(tokens.length < numCols)
                {
                    tokens = (String[]) StsMath.arrayAddArray(tokens, asciiFile.getTokens());
                    if(tokens.length > numCols)
                        StsMath.trimArray(tokens,numCols);
                }
                /*
                for(int j=0; j<tokens.length; j++)
                    row[j] = labelFormat.format(Double.parseDouble(tokens[j]));
                */
                row = tokens;
                fileTable.addRow(row);
                fileTable.setRowType(rowCount, StsTablePanelNew.NOT_EDITABLE);
                rowCount++;

                tokens = asciiFile.getTokens(new String[] {" ",","});
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

    public void setColumnHeadings()
    {
        fileTable.removeAllColumns();
        Object[] colNames = new Object[numCols];
        for(int i=0; i<numCols; i++)
            colNames[i] = "Unknown";

        if(xCol < numCols)
            colNames[xCol] = "X";
        if(yCol < numCols)
            colNames[yCol] = "Y";
        if(tCol < numCols)
            colNames[tCol] = "T";
        if(ilineCol < numCols)
            colNames[ilineCol] = "InLine";
        if(xlineCol < numCols)
            colNames[xlineCol] = "XLine";

        fileTable.addColumns(colNames);
        fileTable.setPreferredSize(new Dimension(colNames.length * 50, 200));
        }
}
