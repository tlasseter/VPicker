package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
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

public class StsSensorDefinePanel extends JPanel implements ListSelectionListener, ActionListener
{
    private StsSensorLoadWizard wizard;
    private StsSensorDefine wizardStep;

    JPanel jPanel1 = new JPanel();
    ButtonGroup volumeGrp = new ButtonGroup();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JButton dirBtn = new JButton();
    JTextField filenameTxt = new JTextField();
    File newFile = null;
    byte timeType = StsSensorKeywordIO.TIME_ONLY;
    int timeIdx = -1;
    int attIndices[] = new int[] {1};
    int nTokens = 0;
    int rowsInFile = 0;
    String attString = "0";
    private StsTablePanel segyTable = new StsTablePanel();
    private JScrollPane jScrollPane1 = new JScrollPane();

    private JFileChooser chooseDirectory = null;
    JLabel jLabel1 = new JLabel();
    JTextField timeColTxt = new JTextField();
    JLabel jLabel2 = new JLabel();
    StsStringFieldBean attributeColTxt = new StsStringFieldBean(this, "attributeIndicesString", "0",true,"Attribute Columns:");
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
    JLabel jLabel4 = new JLabel();
    StsFloatFieldBean xField = new StsFloatFieldBean();
    StsFloatFieldBean yField = new StsFloatFieldBean();
    StsFloatFieldBean zField = new StsFloatFieldBean();
    float xLoc=0.0f, yLoc=0.0f, zLoc=0.0f;
    int xIdx = -1, yIdx = -1, zIdx = -1;

    StsCheckbox dynamicPositionField = new StsCheckbox("Dynamic", "Get Dynamic Position from File?");

    GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsSensorDefinePanel(StsWizard wizard, StsWizardStep wizardStep) {
        this.wizard = (StsSensorLoadWizard) wizard;
        this.wizardStep = (StsSensorDefine) wizardStep;

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
        xField.initialize(this, "XLoc", true, "X:");
        xField.setEditable(true);
        yField.initialize(this, "YLoc", true, "Y:");
        yField.setEditable(true);
        zField.initialize(this, "ZLoc", true, "Z:");
        zField.setEditable(true);
        dynamicPositionField.setSelected(false);
        dynamicPositionField.addActionListener(this);
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
        attributeColTxt.setText("1");
        timeToggle.setMargin(new Insets(2, 2, 2, 2));
        timeToggle.setText("Select");
        jPanel2.setLayout(gridBagLayout2);
        jPanel3.setLayout(gridBagLayout3);
        jPanel3.setBorder(null);
        jPanel2.setBorder(null);
        startTimeTxt.setText("01-01-71 00:00:00.0");
        jLabel3.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel3.setText("Start Date/Time:");
        jLabel4.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel4.setText("Selected File:");

        this.add(jPanel1, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
                                                 , GridBagConstraints.CENTER,
                                                 GridBagConstraints.BOTH,
                                                 new Insets(5, 5, 5, 5), 0, 0));

        jScrollPane1.getViewport().add(segyTable, null);
        jPanel1.add(jScrollPane1,  new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(6, 1, 4, 4), -621, -787));
        jPanel1.add(jPanel2,    new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(1, 1, 0, 4), 3, 3));
        jPanel2.add(filenameTxt,  new GridBagConstraints(1, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 1, 3, 0), 0, 0));
        jPanel2.add(jLabel4,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 21, 3, 0), 0, 0));
        jPanel2.add(dirBtn,  new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 0, 3, 0), 12, 0));
        jPanel2.add(startTimeTxt,   new GridBagConstraints(1, 1, 3, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 0, 2, 0), 0, 0));
        jPanel2.add(jLabel3,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 21, 3, 0), 0, 0));
        jPanel1.add(jPanel3,    new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 1, 0, 4), 10, 0));

        jPanel3.add(attributeColTxt,  new GridBagConstraints(0, 0, 7, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 0, 4, 0), 0, 0));
        jPanel3.add(xField,  new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 4, 0), 0, 0));
        jPanel3.add(yField,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 4, 0), 0, 0));
        jPanel3.add(zField,  new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 4, 0), 0, 0));
        jPanel3.add(dynamicPositionField,  new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 4, 0), 0, 0));

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
            {
                computeTimeIdx();
                if(timeIdx < 0)
                {
                    new StsMessage(wizard.frame, StsMessage.ERROR, "Selected file must have a column labeled TIME");
                    filenameTxt.setText("Invalid File Selected, try again");
                    return;
                }
                fillTable();
            }
        }
        else if(source == dynamicPositionField)
        {
            if(dynamicPositionField.isSelected())
            {
                if (!canBeDynamic())
                {
                    new StsMessage(wizard.frame, StsMessage.ERROR,
                        "Dynamic not possible. X, Y and Depth do not exist in file");
                    dynamicPositionField.setSelected(false);
                    return;
                }
            }
            if(dynamicPositionField.isSelected())
            {
                xField.setEditable(false);
                yField.setEditable(false);
                zField.setEditable(false);
            }
            else
            {
                xField.setEditable(true);
                yField.setEditable(true);
                zField.setEditable(true);
            }
        }
    }

    public boolean isDynamic() { return dynamicPositionField.isSelected(); }
    public void setXLoc(float x) { this.xLoc = x; }
    public void setYLoc(float y) { this.yLoc = y; }
    public void setZLoc(float z) { this.zLoc = z; }
    public float getXLoc() { return xLoc; }
    public float getYLoc() { return yLoc; }
    public float getZLoc() { return zLoc; }
    public int getTimeIndex() { return timeIdx; }
    public void valueChanged(ListSelectionEvent e)
    {
        if (segyTable.getNumberOfRows() > 0)
        {
            if(segyTable.getSelectedColumns().length == 0)
                return;
            attIndices = segyTable.getSelectedColumns();
            attString = null;
            for(int i=0; i<attIndices.length; i++)
            {
                if(attIndices[i] != timeIdx)
                {
                    if(attString != null)
                        attString = attString + "," + (attIndices[i]+1);
                    else
                        attString = new Integer(attIndices[i]+1).toString();
                }
            }
            attributeColTxt.setText(attString);
            if(attString != null)
                getAttributeIndices();
        }
    }

    public void computeTimeIdx()
    {
        BufferedReader bufRdr = null;
        String line = null;
        StringTokenizer stok = null;
        timeIdx = -1;
        try
        {
            bufRdr = new BufferedReader(new FileReader(newFile));
            line = bufRdr.readLine();
            stok = new StringTokenizer(line, ", ;");
            nTokens = stok.countTokens();
            for (int i = 0; i < nTokens; i++)
            {
                if (stok.nextToken().equalsIgnoreCase(StsSensorKeywordIO.TIME))
                {
                    timeIdx = i;
                    break;
                }
            }
        }
        catch(Exception e)
        {
            StsMessageFiles.logMessage("Failed to find time column in file: " + newFile.getName());
        }
    }

    public void determineTimeType()
    {
        BufferedReader bufRdr = null;
        String line = null;
        StringTokenizer stok = null;
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
                    if(time.indexOf("-") < 0)
                    {
                        if(time.indexOf(":") < 0)
                            timeType = StsSensorKeywordIO.ELAPSED_TIME;
                        else
                            timeType = StsSensorKeywordIO.TIME_ONLY;
                        startTimeTxt.setEnabled(true);
                    }
                    else
                    {
                        timeType = StsSensorKeywordIO.TIME_AND_DATE;
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

    public boolean canBeDynamic()
    {
        BufferedReader bufRdr = null;
        String line = null;
        StringTokenizer stok = null;
        int nTokens = 0;
        try
        {
            bufRdr = new BufferedReader(new FileReader(newFile));
            line = bufRdr.readLine();
            stok = new StringTokenizer(line,", ;");
            nTokens = stok.countTokens();
            for(int i=0; i<nTokens; i++)
            {
                String attribute = stok.nextToken();
                if(attribute.equalsIgnoreCase(StsSensorKeywordIO.X))
                    xIdx = i;
                else if(attribute.equalsIgnoreCase(StsSensorKeywordIO.Y))
                    yIdx = i;
                else
                {
                    for(int j=0; j<StsSensorKeywordIO.Z_KEYWORDS.length; j++)
                    {
                        if(attribute.equalsIgnoreCase(StsSensorKeywordIO.Z_KEYWORDS[j]))
                            zIdx = j;
                    }
                }
            }
            bufRdr.close();
            if((xIdx != -1) && (yIdx != -1) && (zIdx != -1))
                return true;
            return false;
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Failed to read file: " + newFile.getName());
            return false;
        }
    }

    public byte getPositionType()
    {
        if(dynamicPositionField.isSelected())
            return StsSensor.DYNAMIC;
        else
        {
            if((xLoc == 0.0f) && (yLoc == 0.0f) && (zLoc == 0.0f))
                return StsSensor.NONE;
            else
                return StsSensor.STATIC;
        }
    }

    public double[] getPosition() { return new double[] {getXLoc(), getYLoc(), getZLoc() }; }

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

        if(getAttributeIndicesString() == null) return null;

        stok = new StringTokenizer(getAttributeIndicesString(),", ;");
        int nTokens = stok.countTokens();
        attIndices = new int[nTokens];
        for(int i=0; i<nTokens; i++)
        {
            token = stok.nextToken();
            attIndices[i] = Integer.parseInt(token) - 1;
        }
        Arrays.sort(attIndices);
        return attIndices;
    }

    public String getAttributeIndicesString()
    {
        return attString;
    }

    public void setAttributeIndicesString(String attStg)
    {
        this.attString = attStg;
    }

    public int getNumberTokens() { return nTokens; }
    public int getNumberValidRows() { return rowsInFile; }
    public String getStartTime()
    {
        Date date = new Date(0L);
        SimpleDateFormat format = new SimpleDateFormat(wizard.getModel().getProject().getTimeDateFormatString());
        String time = format.format(date);
        try
        {
            time = startTimeTxt.getText();
            time = wizard.cleanTimeString(time);
            if(time == null)
            {
                new StsMessage(wizard.frame, StsMessage.ERROR,
                    "Invalid date & time value (" + format.format(date) + ") in time series file selection step.\n" +
                    "\n   Solution: Return to step and re-enter valid start time.");
                return null;
            }
            date = format.parse(time);
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Failed to create date, setting to 01-01-71 00:00:00.0");
            return "01-01-71 00:00:00.0";
        }
        return time;
    }

    public void verifyXYZSelection()
    {
        int[] indices = null;
        getAttributeIndices();
        if(StsMath.binarySearch(attIndices, xIdx) < 0)
        {
            StsMessageFiles.infoMessage("X Column required for dynamic sensors, added to user list.");
            indices = new int[attIndices.length + 1];
            System.arraycopy(attIndices, 0, indices, 0, attIndices.length);
            indices[attIndices.length] = xIdx;
            attIndices = indices;
            Arrays.sort(attIndices);
        }
        if(StsMath.binarySearch(attIndices, yIdx) < 0)
        {
            StsMessageFiles.infoMessage("Y Column required for dynamic sensors, added to user list.");
            indices = new int[attIndices.length + 1];
            System.arraycopy(attIndices, 0, indices, 0, attIndices.length);
            indices[attIndices.length] = yIdx;
            attIndices = indices;
            Arrays.sort(attIndices);
        }
        if(StsMath.binarySearch(attIndices, zIdx) < 0)
        {
            StsMessageFiles.infoMessage("Z Column required for dynamic sensors, added to user list.");
            indices = new int[attIndices.length + 1];
            System.arraycopy(attIndices, 0, indices, 0, attIndices.length);
            indices[attIndices.length] = zIdx;
            attIndices = indices;
            Arrays.sort(attIndices);
        }

        // Build Attribute String
        String attString = new String();
        for(int i=0; i<attIndices.length; i++)
        {
            if(i < attIndices.length-1)
                attString = attString + (attIndices[i]+1) + ",";
            else
                attString = attString + (attIndices[i]+1);
        }
        attributeColTxt.setText(attString);
    }

    public boolean verifyUserAttributes()
    {
        // Are there any available columns to pick beyond those required?
        if(isDynamic())
        {
            if(nTokens == 4)
                return true;
        }
        for(int i=0; i<attIndices.length; i++)
        {
            if((attIndices[i] != zIdx) && (attIndices[i] != xIdx) && (attIndices[i] != yIdx) && (attIndices[i] != timeIdx))
                return true;
        }
        return false;
    }
}
