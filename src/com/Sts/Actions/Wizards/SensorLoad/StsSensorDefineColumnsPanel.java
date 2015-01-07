package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Import.StsKeywordIO;
import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.DBTypes.StsSensor;
import com.Sts.IO.StsAsciiFile;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Beans.StsListFieldBean;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsTablePanel;
import com.Sts.UI.StsTablePanelNew;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsParameters;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSensorDefineColumnsPanel extends StsJPanel implements ActionListener, ChangeListener
{
    private StsTablePanel fileTable = new StsTablePanel();
	private JScrollPane tableScrollPane1 = new JScrollPane();
    private JLabel numColLabel = new JLabel("# Cols:");

    StsListFieldBean fileListBean;
    private Object[] selectedFiles = null;

    private StsGroupBox colDefPanel = new StsGroupBox("Column Assign");
    private StsSensorFile currentFile = null;
    private StsSensorLoadWizard wizard;
    private StsSensorDefineColumns wizardStep;

    int numSkipped = 0;

    JSpinner numColSpin = new JSpinner();
    private SpinnerModel numColSpinModel = null;
    JSpinner xSpin = new JSpinner();
    private SpinnerModel xSpinModel = null;
    JSpinner ySpin = new JSpinner();
    private SpinnerModel ySpinModel = null;
    JSpinner zSpin = new JSpinner();
    private SpinnerModel zSpinModel = null;
    JSpinner timeSpin = new JSpinner();
    private SpinnerModel timeSpinModel = null;
    JSpinner dateSpin = new JSpinner();
    private SpinnerModel dateSpinModel = null;
    JLabel jLabel3 = new JLabel();
    JLabel jLabel7 = new JLabel();
    JLabel jLabel8 = new JLabel();
    JLabel jLabel9 = new JLabel();
    JLabel jLabel10 = new JLabel();
    JLabel jLabel11 = new JLabel();
    JLabel jLabel12 = new JLabel();
    JCheckBox dateTimeChk = new JCheckBox("Date & Time are in the same column");
    boolean dateTime = true;

    JPanel jPanel2 = new JPanel();
    JButton okBtn = new JButton();
    JButton cancelBtn = new JButton();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    public StsSensorDefineColumnsPanel(StsWizard wizard, StsSensorDefineColumns wizardStep)
    {
        this.wizard = (StsSensorLoadWizard)wizard;
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

        fileListBean = new StsListFieldBean();
        JScrollPane fileScrollPane = new JScrollPane();
        fileScrollPane.getViewport().add(fileListBean, null);

        StsJPanel fileBox = new StsGroupBox("Select File");
        fileBox.setPreferredSize(300, 150);
        fileBox.gbc.fill = gbc.BOTH;
        fileBox.gbc.weighty = 1.0;
        fileBox.addEndRow(fileScrollPane);
        gbc.fill = gbc.BOTH;
        gbc.weighty = 0.2f;
        add(fileBox);

        numColSpinModel = new SpinnerNumberModel(5, 1, 200, 1);
        xSpinModel = new SpinnerNumberModel(1, 1, 200, 1);
        ySpinModel = new SpinnerNumberModel(2, 1, 200, 1);
        zSpinModel = new SpinnerNumberModel(3, 1, 200, 1);
        timeSpinModel = new SpinnerNumberModel(4, 1, 200, 1);
        dateSpinModel = new SpinnerNumberModel(5, 1, 200, 1);

        fileTable.setTitle("Sensor File:");
		tableScrollPane1.setAutoscrolls(true);
		tableScrollPane1.getViewport().add(fileTable, null);
		tableScrollPane1.getViewport().setPreferredSize(new Dimension(550, 100));

        jLabel3.setFont(new Font("Serif", 1, 11));
        jLabel3.setHorizontalAlignment(SwingConstants.CENTER);
        jLabel3.setText("Field");
        jLabel7.setFont(new Font("Serif", 1, 11));
        jLabel7.setHorizontalAlignment(SwingConstants.CENTER);
        jLabel7.setText("Column");
        jLabel8.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel8.setText("X:");
        jLabel9.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel9.setText("Y:");
        jLabel10.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel10.setText("Z:");
        jLabel11.setText("Time:");
        jLabel12.setText("Date:");
        numColSpin.addChangeListener(this);
        numColSpin.setModel(numColSpinModel);
        xSpin.addChangeListener(this);
        xSpin.setModel(xSpinModel);
        ySpin.addChangeListener(this);
        ySpin.setModel(ySpinModel);
        zSpin.addChangeListener(this);
        zSpin.setModel(zSpinModel);
        timeSpin.addChangeListener(this);
        timeSpin.setModel(timeSpinModel);
        dateSpin.addChangeListener(this);
        dateSpin.setModel(dateSpinModel);
        dateTimeChk.addActionListener(this);
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout3);
        okBtn.setText("Ok");
        okBtn.addActionListener(this);
        cancelBtn.setText("Cancel");
        cancelBtn.addActionListener(this);

        colDefPanel.addToRow(numColLabel);
        colDefPanel.addEndRow(numColSpin);
        colDefPanel.gbc.gridwidth = 2;
        colDefPanel.addEndRow(dateTimeChk);
        colDefPanel.gbc.gridwidth = 1;
        colDefPanel.addToRow(jLabel8);
        colDefPanel.addEndRow(xSpin);
        colDefPanel.addToRow(jLabel9);
        colDefPanel.addEndRow(ySpin);
        colDefPanel.addToRow(jLabel10);
        colDefPanel.addEndRow(zSpin);
        colDefPanel.addToRow(jLabel11);
        colDefPanel.addEndRow(timeSpin);
        colDefPanel.addToRow(jLabel12);
        colDefPanel.addEndRow(dateSpin);
        dateSpin.setEnabled(false);

        StsGroupBox colViewPanel = new StsGroupBox("Selected File Column Definition");
        colViewPanel.gbc.weightx = 0.0;
        colViewPanel.gbc.fill = fileBox.gbc.VERTICAL;
        colViewPanel.gbc.gridwidth = 1;
        colViewPanel.addToRow(colDefPanel);
        colViewPanel.gbc.weightx = 1.0;
        colViewPanel.gbc.fill = fileBox.gbc.BOTH;
        colViewPanel.addEndRow(tableScrollPane1);
        gbc.weighty = 0.8f;
        add(colViewPanel);
    }

    public void initialize()
    {
        StsSensorFile[] sensorFiles = wizard.getSensorFiles();
        for(int i=0; i<sensorFiles.length; i++)
        {
            // No Header in File.
            if(sensorFiles[i].getNumCols() == 0)
                sensorFiles[i].setNumCols(2); // Minimum number of required columns is 2
        }
        numSkipped = sensorFiles[0].getNumHeaderRows();
        currentFile = sensorFiles[0];
        fileListBean.initialize(this,"files",null,sensorFiles);
        fileListBean.setSelectedIndex(0);
        selectedFiles = new Object[] {currentFile};

        numColSpin.setValue(new Integer(currentFile.getNumCols()));

        xSpin.setValue(new Integer(currentFile.colOrder[currentFile.X]+1));
        ySpin.setValue(new Integer(currentFile.colOrder[currentFile.Y]+1));
        zSpin.setValue(new Integer(currentFile.colOrder[currentFile.Z]+1));
        timeSpin.setValue(new Integer(currentFile.colOrder[currentFile.TIME]+1));
        dateSpin.setValue(new Integer(currentFile.colOrder[currentFile.DATE]+1));

        dateTimeChk.setSelected(dateTime);
        wizard.rebuild();
        updateFileView();

        fileListBean.setSelectedAll();
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == dateTimeChk)
        {
        	dateTime = dateTimeChk.isSelected();
            if(dateTime)
                dateSpin.setEnabled(true);
            else
                dateSpin.setEnabled(false);
        }
    }

    public void stateChanged(ChangeEvent e)
    {
        Object source = e.getSource();
        if(source == xSpin)
            setLocations(StsSensorFile.X, Integer.valueOf(xSpin.getValue().toString()).intValue() - 1);
        else if(source == ySpin)
            setLocations(StsSensorFile.Y, Integer.valueOf(ySpin.getValue().toString()).intValue() - 1);
        else if(source == zSpin)
            setLocations(StsSensorFile.Z, Integer.valueOf(zSpin.getValue().toString()).intValue() - 1);
        else if(source == timeSpin)
            setLocations(StsSensorFile.TIME, Integer.valueOf(timeSpin.getValue().toString()).intValue() - 1);
        else if(source == dateSpin)
            setLocations(StsSensorFile.DATE, Integer.valueOf(dateSpin.getValue().toString()).intValue() - 1);
        else if(source == numColSpin)
        {
            Object[] files = fileListBean.getSelectedObjects();
            for(int i=0; i<files.length; i++)
                ((StsSensorFile)files[i]).setNumCols(Integer.valueOf(numColSpin.getValue().toString()).intValue());
        }
        if(currentFile.getNumCols() > 0)
            updateFileView();
    }

    public void setLocations(byte type, int col)
    {
        if(selectedFiles == null) return;        
        for(int i=0; i<selectedFiles.length; i++)
           ((StsSensorFile)selectedFiles[i]).setColLocation(type,col);
    }
    public void setFiles(Object file)
    {
        selectedFiles = fileListBean.getSelectedObjects();
        currentFile = (StsSensorFile)file;

        numColSpin.setValue(new Integer(currentFile.getNumCols()));

        xSpin.setValue(new Integer(currentFile.colOrder[currentFile.X]+1));
        ySpin.setValue(new Integer(currentFile.colOrder[currentFile.Y]+1));
        zSpin.setValue(new Integer(currentFile.colOrder[currentFile.Z]+1));
        timeSpin.setValue(new Integer(currentFile.colOrder[currentFile.TIME]+1));
        dateSpin.setValue(new Integer(currentFile.colOrder[currentFile.DATE]+1));

        if(currentFile.positionType == StsSensor.STATIC)
        {
            xSpin.setEnabled(false);
            ySpin.setEnabled(false);
            zSpin.setEnabled(false);
        }
        else
        {
            xSpin.setEnabled(true);
            ySpin.setEnabled(true);
            zSpin.setEnabled(true);
        }
        updateFileView();
    }

    public Object getFiles() { return currentFile; }

    public void skipHeader(StsAsciiFile asciiFile, boolean determineOrder)
    {
        String line = null;
        try
        {
            for (int i = 0; i < currentFile.numHeaderRows; i++)
                line = asciiFile.readLine();
        }
        catch(Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING,  e.getMessage());
            return;
        }
    }

    public boolean dateTime() { return dateTime; }

    public boolean updateFileView()
    {
        String[] tokens;

        StsAsciiFile asciiFile = new StsAsciiFile(currentFile.file);
        if(!asciiFile.openReadWithErrorMessage()) return false;

        try
        {
            skipHeader(asciiFile, false);
            tokens = asciiFile.getTokens(currentFile.delimiters);
            int nLines = 0;
            Object[] row = new Object[currentFile.getNumCols()];

            fileTable.removeAllRows();
            setColumnHeadings();

            int rowCount = 0;
            while(tokens != null)
            {
                if(++nLines > 100)
                    break;

                while(tokens.length < currentFile.getNumCols())
                {
                    tokens = (String[]) StsMath.arrayAddArray(tokens, asciiFile.getTokens());
                    if(tokens.length > currentFile.getNumCols())
                        StsMath.trimArray(tokens,currentFile.getNumCols());
                }
                /*
                for(int j=0; j<tokens.length; j++)
                    row[j] = labelFormat.format(Double.parseDouble(tokens[j]));
                */
                row = tokens;
                fileTable.addRow(row);
                fileTable.setRowType(rowCount, StsTablePanelNew.NOT_EDITABLE);
                rowCount++;

                tokens = asciiFile.getTokens(currentFile.delimiters);
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

    public void setColumnHeadings()
    {
        fileTable.removeAllColumns();
        Object[] colNames = new Object[currentFile.getNumCols()];
        if(currentFile.getNumCols() <= 0) return;
        int cnt = 0;
        for(int i=0; i<currentFile.getNumCols(); i++)
        {
            if(currentFile.curveNames != null)
            {
                if(i == currentFile.getColLocation(currentFile.TIME))
                {
                    colNames[i] = "TIME";
                    cnt++;
                }
                else if((currentFile.getColLocation(currentFile.DATE) != -1) && (i == currentFile.getColLocation(currentFile.DATE)))
                {
                    colNames[i] = "DATE";
                    cnt++;
                }
                else
                    colNames[i] = currentFile.curveNames[i-cnt];
            }
            else
            {
                  colNames[i] = "Unknown";
            }
        }

        if(currentFile.getColLocation(currentFile.X) < currentFile.getNumCols() && (currentFile.getColLocation(currentFile.X) != -1))
            colNames[currentFile.getColLocation(currentFile.X)] = "X";
        if(currentFile.getColLocation(currentFile.Y) < currentFile.getNumCols() && (currentFile.getColLocation(currentFile.Y) != -1))
            colNames[currentFile.getColLocation(currentFile.Y)] = "Y";
        if(currentFile.getColLocation(currentFile.Z) < currentFile.getNumCols() && (currentFile.getColLocation(currentFile.Z) != -1))
            colNames[currentFile.getColLocation(currentFile.Z)] = "Z";
        if(currentFile.getColLocation(currentFile.TIME) < currentFile.getNumCols() && (currentFile.getColLocation(currentFile.TIME) != -1))
            colNames[currentFile.getColLocation(currentFile.TIME)] = "Time";
        if((currentFile.getColLocation(currentFile.DATE) < currentFile.getNumCols()) && (currentFile.getColLocation(currentFile.DATE) != -1))
        {
            colNames[currentFile.getColLocation(currentFile.DATE)] = "Date";
            dateSpin.setEnabled(true);
        }
        else
            dateSpin.setEnabled(false);


        fileTable.addColumns(colNames);
//        fileTable.setPreferredSize(new Dimension(colNames.length * 50, 200));
        }
}