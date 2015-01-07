package com.Sts.Actions.Wizards.TimeSeries;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
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

public class StsTimeSeriesAsciiPanel extends JPanel implements ListSelectionListener, ActionListener
{
    private StsTimeSeriesWizard wizard;
    private StsTimeSeriesAscii wizardStep;

    JPanel jPanel1 = new JPanel();
    ButtonGroup volumeGrp = new ButtonGroup();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JButton dirBtn = new JButton();
    JTextField filenameTxt = new JTextField();
    File newFile = null;
    byte timeType = wizard.ARTIFICIAL_TIME;
    int timeIdx = 0;
    int attIndices[] = new int[] {1};
    int nTokens = 0;
    int rowsInFile = 0;

    private StsTablePanel segyTable = new StsTablePanel();
    private JScrollPane jScrollPane1 = new JScrollPane();

    private JFileChooser chooseDirectory = null;
    JLabel jLabel1 = new JLabel();
    JTextField timeColTxt = new JTextField();
    JLabel jLabel2 = new JLabel();
    JTextField attributeColTxt = new JTextField();
    JToggleButton timeToggle = new JToggleButton();
    JToggleButton attributeToggle = new JToggleButton();
    ButtonGroup timeAttGroup = new ButtonGroup();
    JPanel jPanel2 = new JPanel();
    JPanel jPanel3 = new JPanel();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    GridBagLayout gridBagLayout4 = new GridBagLayout();
    ButtonGroup elapseClockGrp = new ButtonGroup();
    JTextField startTimeTxt = new JTextField();
    JLabel jLabel3 = new JLabel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsTimeSeriesAsciiPanel(StsWizard wizard, StsWizardStep wizardStep) {
        this.wizard = (StsTimeSeriesWizard) wizard;
        this.wizardStep = (StsTimeSeriesAscii) wizardStep;

        try {
            jbInit();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        setPreferredSize(new java.awt.Dimension(400, 400));
        wizard.rebuild();
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout1);
        jPanel1.setLayout(gridBagLayout4);

        segyTable.setTitle("File Contents:");
        segyTable.setFont(new Font("Dialog", 3, 12));
        segyTable.setColumnSelectable(true);
        segyTable.setRowsSelectable(false);
        segyTable.addListSelectionListener(this);
        jScrollPane1.setAutoscrolls(true);
        jScrollPane1.setPreferredSize(new Dimension(1000, 1000));

        dirBtn.setMargin(new Insets(2, 2, 2, 2));
        dirBtn.setText("Select...");
        dirBtn.addActionListener(this);

        filenameTxt.setText("Ascii File");
        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setText("Time Column:");
        timeColTxt.setText("0");
        jLabel2.setText("Attribute Columns:");
        jLabel2.setHorizontalAlignment(SwingConstants.RIGHT);
        attributeColTxt.setText("1");
//        timeToggle.setMargin(new Insets(2, 2, 2, 2));
//        timeToggle.setText("Select");
//        attributeToggle.setText("Select");
//        attributeToggle.setMargin(new Insets(2, 2, 2, 2));
        jPanel2.setLayout(gridBagLayout2);
        jPanel3.setLayout(gridBagLayout3);
        jPanel3.setBorder(null);
        jPanel2.setBorder(null);
        startTimeTxt.setText("12:00:00");
        jLabel3.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel3.setText("Start Time:");
        this.add(jPanel1, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
                                                 , GridBagConstraints.CENTER,
                                                 GridBagConstraints.BOTH,
                                                 new Insets(5, 5, 5, 5), 0, 0));

        jScrollPane1.getViewport().add(segyTable, null);
        jPanel1.add(jScrollPane1,  new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(6, 1, 4, 4), -621, -787));
        jPanel1.add(jPanel2,    new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(1, 1, 0, 4), 3, 3));
        jPanel2.add(filenameTxt,  new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 1, 3, 0), 111, 3));
        jPanel2.add(dirBtn,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 0, 3, 0), 12, 0));
        jPanel2.add(startTimeTxt,   new GridBagConstraints(3, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 2, 0), 28, 1));
        jPanel2.add(jLabel3,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 21, 3, 0), 5, 7));
        jPanel1.add(jPanel3,    new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 1, 0, 4), 10, 0));
        jPanel3.add(jLabel1,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 1, 4, 0), 5, 7));
        jPanel3.add(timeColTxt,  new GridBagConstraints(1, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 0, 4, 0), 16, 1));
//        jPanel3.add(timeToggle,  new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0
//            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(2, 0, 4, 0), 16, 0));
        jPanel3.add(jLabel2,  new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 4, 0), 9, 7));
        jPanel3.add(attributeColTxt,  new GridBagConstraints(4, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 0, 4, 0), 61, 1));
//        jPanel3.add(attributeToggle,  new GridBagConstraints(5, 0, 1, 1, 1.0, 1.0
//            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(2, 0, 4, 1), 7, 0));

        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        segyTable.addListSelectionListener(this);
        timeAttGroup.add(timeToggle);
        timeAttGroup.add(attributeToggle);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == dirBtn)
        {
            directoryBrowse();
            if(newFile != null)
                fillTable();
        }
    }

    public void valueChanged(ListSelectionEvent e)
    {

        if (segyTable.getNumberOfRows() > 0)
        {
            if(segyTable.getSelectedColumns().length == 0)
                return;
            attIndices = segyTable.getSelectedColumns();
            if(timeToggle.isSelected())
            {
                timeColTxt.setText(new Integer(attIndices[0]).toString());
                timeIdx = attIndices[0];
                determineTimeType();
            }
            else
            {
                String attString = null;
                for(int i=0; i<attIndices.length; i++)
                {
                    if(i != 0)
                        attString = attString + "," + attIndices[i];
                    else
                        attString = new Integer(attIndices[0]).toString();
                }
                attributeColTxt.setText(attString);
            }
        }
    }

    public void determineTimeType()
    {
        BufferedReader bufRdr = null;
        String line = null;
        StringTokenizer stok = null;
        timeIdx = getTimeIndex();
        int nTokens = 0;
        try
        {
            bufRdr = new BufferedReader(new FileReader(newFile));
            line = bufRdr.readLine();
            while((line = bufRdr.readLine()) != null)
            {
                stok = new StringTokenizer(line,", ;");
                nTokens = stok.countTokens();
                if(nTokens < this.nTokens)
                    continue;
                if(nTokens < timeIdx+1)
                    continue;
                else
                    break;
            }
            String time = null;
            for(int i=0; i<nTokens; i++)
            {
                if(i == timeIdx)
                {
                    time = stok.nextToken();
                    if(time.indexOf(":") < 0)
                    {
                        timeType = wizard.ARTIFICIAL_TIME;
                        startTimeTxt.setEnabled(true);
                    }
                    else
                    {
                        timeType = wizard.ACTUAL_TIME;
                        startTimeTxt.setEnabled(false);
                    }
                }
                else
                    stok.nextToken();
            }
            bufRdr.close();
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Failed to read file: " + newFile.getName());
            return;
        }
    }

    private void fillTable()
    {
        BufferedReader bufRdr = null;
        String line = null;
        StringTokenizer stok = null;
        boolean firstRow = true;
        int rowCount = 0;
        Object[] row = null;

        segyTable.removeAllColumns();
        segyTable.removeAllRows();
        try
        {
            bufRdr = new BufferedReader(new FileReader(newFile));
            while(((line = bufRdr.readLine()) != null) || (rowCount == 25))
            {
                stok = new StringTokenizer(line,", ;");
                nTokens = stok.countTokens();
                if((firstRow) && (nTokens < 2))
                {
                    new StsMessage(wizard.frame, StsMessage.ERROR, "Comma delimited file must have at least 2 columns.");
                    return;
                }

                else if(nTokens < 2)
                {
                    rowsInFile--;
                    continue;
                }

                // Create column labels if exist in file
                if(firstRow)
                {
                    firstRow = false;
                    segyTable.setPreferredSize(new Dimension(nTokens * 50, 500));
                    segyTable.setSize(new Dimension(nTokens * 50, 500));
                    for(int i=0; i<nTokens; i++)
                    {
                        String token = stok.nextToken();
                        if(StsStringUtils.isNumeric(token))
                            token = new String("C" + i);

                        if(isAttribute(i))
                            segyTable.addAccentedColumn(token);
                        else
                            segyTable.addColumn(token);
                    }
                    row = new Object[nTokens];
                    rowsInFile--;
                    continue;
                }
                for(int i=0; i<nTokens; i++)
                {
                    row[i] = stok.nextElement();
                }
                segyTable.addRow(row);
                rowCount++;
                rowsInFile++;
            }
            while((line = bufRdr.readLine()) != null)
            {
                stok = new StringTokenizer(line,", ;");
                if(nTokens == stok.countTokens())
                    rowsInFile++;
            }
            bufRdr.close();
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Failed to read file: " + newFile.getName());
            return;
        }
    }

    private boolean isAttribute(int idx)
    {
        if(idx == timeIdx)
            return true;
        if(attIndices == null)
            return false;
        for(int i=0; i<attIndices.length; i++)
            if(idx == attIndices[i])
                return true;
        return false;
    }

    public void directoryBrowse()
    {
        String currentDirectory = wizard.getModel().getProject().getRootDirString();
        chooseDirectory = new JFileChooser(currentDirectory);
        chooseDirectory.setFileSelectionMode(JFileChooser.FILES_ONLY);
        chooseDirectory.setDialogTitle("Select or Enter Desired File and Press Open");
        chooseDirectory.setApproveButtonText("Open File");

        chooseDirectory.showOpenDialog(null);
        newFile = chooseDirectory.getSelectedFile();
        if(newFile != null)
            filenameTxt.setText(newFile.getName());
    }
    public File getAsciiFile() { return newFile; }
    public byte getAsciiTimeType() { return timeType; }
    public int[] getAttributeIndices()
    {
        String token = null;
        StringTokenizer stok = null;
        stok = new StringTokenizer(attributeColTxt.getText(),", ;");
        int nTokens = stok.countTokens();
        int[] atts = new int[nTokens];
        for(int i=0; i<nTokens; i++)
        {
            token = stok.nextToken();
            atts[i] = Integer.parseInt(token);
        }
        return atts;
    }
    public int getTimeIndex()
    {
        String token = timeColTxt.getText();
        return Integer.parseInt(token);
    }
    public int getNumberTokens() { return nTokens; }
    public int getNumberValidRows() { return rowsInFile; }
    public long getStartTime()
    {
        Date date = new Date(0L);
        try
        {
            SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss.S");
            String time = startTimeTxt.getText();
            time = wizard.cleanTimeString(time);
            if(time == null)
            {
                new StsMessage(wizard.frame, StsMessage.ERROR,
                    "Invalid time value (##:##:##) in time series file selection step.\n" +
                    "\n   Solution: Return to step and re-enter valid start time.");
                return 0L;
            }
            date = format.parse(time);
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Failed to create date, setting to 0L");
            return 0L;
        }
        return date.getTime();
    }
}
