package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsWellSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
    private StsFileTransferPanelNew selectionPanel;

    String[] formats = new String[] {"ut","las","txt"};
    String groups = "well-dev";
    StsFilenameFilter filenameFilter  = new WellFilenameFilter(groups);

    private String currentDirectory = null;
    private String outputDirectory = null;
    private StsAbstractFile currentSelectedFile = null;

    private boolean archiveIt = false;
    private boolean reloadIt = false;
    private boolean hasReloadButton = true;
    private boolean hasOverrideButton = true;

    StsComboBoxFieldBean hUnitsBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean vUnitsBean = new StsComboBoxFieldBean();
    StsFloatFieldBean datumShiftBean = new StsFloatFieldBean();
    float datumShift = 0.0f;
    StsCheckbox applyToTdBtn = new StsCheckbox("Time-Depth", "Apply the correction to the time-depth data.");
    StsCheckbox applyToRefBtn = new StsCheckbox("Markers", "Apply the correction to the marker data.");
    StsCheckbox applyToLogBtn = new StsCheckbox("Logs", "Apply the correction to the log data.");
    StsCheckbox applyToDevBtn = new StsCheckbox("Path", "Apply the correction to the well path data.");
    StsCheckbox[] applyToBtns = new StsCheckbox[] { applyToDevBtn, applyToLogBtn,  applyToTdBtn, applyToRefBtn };

    private StsWellWizard wizard;
    private StsWellSelect wizardStep;
    private StsFilenameFilter userFilter = null;

    private byte vUnits, hUnits;

    public static int S2S = 0;
    public static int LAS = 1;
    public static int UT = 2;

    ButtonGroup formatGroup = new ButtonGroup();
    StsJPanel unitFormatBox = new StsGroupBox("Units and Format");
    StsJPanel datumShiftBox = new StsGroupBox("Depth Correction");

    private StsButton viewButton = null;

    public StsWellSelectPanel(StsWizard wizard, StsWizardStep wizardStep, boolean reloadAscii)
    {
        this.wizard = (StsWellWizard)wizard;
        this.wizardStep = (StsWellSelect)wizardStep;
        hasReloadButton = reloadAscii;

        try
        {
            if (wizard.getModel() != null)
            {
                StsProject project = wizard.getModel().getProject();
                currentDirectory = project.getRootDirString();
                outputDirectory = project.getDataFullDirString();
                String plugIn = wizard.getModel().getWorkflowPlugIn().name;
                if(plugIn.equals("StsVspWorkflow"))
                    archiveIt = true;
            }
            else
            {
                currentDirectory = System.getProperty("user.dirNo"); // standalone testing
                outputDirectory = currentDirectory;
            }
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void constructBeans()
    {
        wizard.constructAvailableFiles();
        //selectionPanel = new StsFileTransferPanelNew(currentDirectory, filenameFilter, this);
		//viewButton = new StsButton("View Selected File", "Review selected file contents", this, "showFile");
        if(wizard.selectFileType.getType() == wizard.GEOGRAPHIX_WELLS)
        {
            filenameFilter = new StsFilenameEndingFilter(new String[] {"wba","wls"});
            selectionPanel.setFilenameFilter(filenameFilter);
        }
        else
            selectionPanel.setFilenameFilter(filenameFilter);
		hUnitsBean.initialize(this, "horzUnitsString", "Horizontal Units:  ", StsParameters.DIST_STRINGS);
		vUnitsBean.initialize(this, "vertUnitsString", "Vertical Units:  ", StsParameters.DIST_STRINGS);
		datumShiftBean.initialize(this, "datumShift", true, "Correction:  ");
        datumShiftBean.setToolTipText("Correction is added to the MDEPTH or DEPTH values");
		datumShiftBean.setColumns(6);
	}

    public void initialize()
    {
        hUnits = wizard.hUnits;
        vUnits = wizard.vUnits;

        constructBeans();

        hUnitsBean.setSelectedItem(StsParameters.DIST_STRINGS[hUnits]);
		vUnitsBean.setSelectedItem(StsParameters.DIST_STRINGS[vUnits]);
        wizard.rebuild();
     }

    void jbInit() throws Exception
    {
        selectionPanel = new StsFileTransferPanelNew(currentDirectory, filenameFilter, this);
		viewButton = new StsButton("View Selected File", "Review selected file contents", this, "showFile");

        gbc.fill = gbc.BOTH;
        gbc.anchor = gbc.NORTH;
        gbc.weighty = 1.0;
        addEndRow(selectionPanel);

        gbc.weighty = 0.0;
        gbc.fill = gbc.NONE;
        gbc.anchor = gbc.CENTER;
        addEndRow(viewButton);

        gbc.anchor = gbc.NORTH;
        gbc.fill = gbc.HORIZONTAL;
        unitFormatBox.addToRow(vUnitsBean);
        unitFormatBox.addEndRow(hUnitsBean);
        addEndRow(unitFormatBox);

        datumShiftBox.addToRow(datumShiftBean);
        datumShiftBox.addToRow(applyToDevBtn);
        datumShiftBox.addToRow(applyToLogBtn);
        datumShiftBox.addToRow(applyToRefBtn);
        datumShiftBox.addEndRow(applyToTdBtn);
        addEndRow(datumShiftBox);

        //gbc.fill = gbc.BOTH;
        //gbc.weighty = 1.0;
        //addEndRow(selectedWellInfoPanel);
    }

    public void showFile()
    {
        String line;
        if(currentSelectedFile == null) return;
        StsAsciiFile selectedFile = new StsAsciiFile(currentSelectedFile);
        if(selectedFile == null)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "No file selected in list, must be in selected (left) list");
            return;
        }
        if(!selectedFile.openReadWithErrorMessage()) return;

        StsFileViewDialog dialog = new StsFileViewDialog(wizard.frame,"Well File View", false);
        dialog.setVisible(true);
        dialog.setViewTitle("File - " + currentSelectedFile.getFilename());
        try {
            while (true)
            {
                line = selectedFile.readLine();
                if (line == null)
                    break;
                dialog.appendLine(line);
            }
        }
        catch (Exception ex)
        {
            StsMessageFiles.infoMessage("Unable to view selected file.");
        }
        finally
        {
            selectedFile.close();
        }
    }
    public boolean hasDirectorySelection() { return true;  }
    
    public void availableFileSelected(StsAbstractFile selectedFile)
    {
    }
    
    public boolean hasReloadButton()
    {
        return hasReloadButton;
    }

    public boolean hasOverrideButton()
    {
        return hasOverrideButton;
    }

    public void setOverrideFilter(boolean override)
    {
        StsWellImport.setOverrideFilter(override);
        // Override dialog if true
        if(override)
        {
            // Show dialog to add a well
            OverrideFilterDialog wDialog = new OverrideFilterDialog(wizard.frame);
            wDialog.show();
            userFilter = wDialog.getFilter();
            if(userFilter != null)
                selectionPanel.setFilenameFilter(userFilter);
        }
        wizard.constructAvailableFiles();
        selectionPanel.setAvailableFiles();
    }

    public StsFilenameFilter getUserFilter() { return userFilter; }
    public boolean getOverrideFilter()
    {
        return StsWellImport.getOverrideFilter();
    }

    public void setReload(boolean reload)
    {
        StsWellImport.setReloadAscii(reload);
        wizard.constructAvailableFiles();
        selectionPanel.setAvailableFiles();
    }

    public boolean getReload()
    {
        return StsWellImport.getReloadAscii();
    }
    public boolean hasArchiveItButton()
    {
        return true;
    }

    public void setArchiveIt(boolean archive)
    {
        archiveIt = archive;
    }

    public boolean getArchiveIt()
    {
        return archiveIt;
    }
    public void fileSelected(StsAbstractFile selectedFile)
    {
        if(selectedFile == null)
        {
            currentSelectedFile = null;
        }
        if (selectedFile != currentSelectedFile)
        {
            currentSelectedFile = selectedFile;
        }
    }

    public String getSelectedDirectory() { return selectionPanel.getCurrentDirectory(); }
    public void addFiles(StsAbstractFile[] files)
    {
        for (int n = 0; n < files.length; n++)
        {
            if(!wizard.addFile(files[n], outputDirectory))
            {
                selectionPanel.removeFile(files[n]);
                break;
            }
        }
        fileSelected(files[0]);
    }

    public void addAllFiles(StsFile[] files)
    {
        for (int n = 0; n < files.length; n++)
            wizard.addFile(files[n], outputDirectory);
        fileSelected(files[files.length-1]);
    }

    public void removeFiles(StsAbstractFile[] files)
    {
        for (int n = 0; n < files.length; n++)
            wizard.removeFile(files[n]);
    }

    public void removeAllFiles()
    {
        wizard.removeFiles();
    }

    public void resetSelection()
    {
        selectionPanel.removeAllFiles();
    }

    public String[] getSelectedWells()
    {
        StsAbstractFile[] files = wizard.getSelectedFiles();
        String[] fileStems = new String[files.length];
        for(int i=0; i<fileStems.length; i++)
        {
            if(userFilter != null)
                fileStems[i] = userFilter.getFilenameStemName(files[i].getFilename());
            else
            {
                fileStems[i] = filenameFilter.getFilenameStemName(files[i].getFilename());
            }
        }
        return fileStems;
    }

    public float getDatumShift() { return datumShift; }
    public void setDatumShift(float shift) { datumShift = shift; }
    public boolean[] getApplyDatumShift()
    {
        boolean[] apply = new boolean[applyToBtns.length];
        for(int i=0; i<applyToBtns.length; i++)
            apply[i] = applyToBtns[i].isSelected();
        return apply;
    }
    public String getHorzUnitsString() { return StsParameters.DIST_STRINGS[hUnits]; }
    public String getVertUnitsString() { return StsParameters.DIST_STRINGS[vUnits]; }
    public byte getHorzUnits() { return hUnits; }
    public byte getVertUnits() { return vUnits; }
    public void setHorzUnitsString(String unitString)
    {
        hUnits = StsParameters.getDistanceUnitsFromString(unitString);
        wizard.hUnits = hUnits;
    }

    public void setVertUnitsString(String unitString)
    {
        vUnits = StsParameters.getDistanceUnitsFromString(unitString);
        wizard.vUnits = vUnits;
    }

    class WellFilenameFilter extends StsFilenameFilter
    {
        String group = null;
        public WellFilenameFilter(String group)  { this.group = group; }

        public boolean accept(File dir, String name)
        {
            if(StsWellImport.getCurrentDirectory() == null)
            {
                StsWellImport.setCurrentDirectory(dir.getAbsolutePath());
                wizard.constructAvailableFiles();
            }
            else
            {
                if(!StsWellImport.getCurrentDirectory().equals(dir.getAbsolutePath()))
                {
                    StsWellImport.setCurrentDirectory(dir.getAbsolutePath());
                    wizard.constructAvailableFiles();
                }
            }
            StsKeywordIO.parseBinaryFilename(name);
            return StsWellImport.getWellFilenameSet(StsKeywordIO.name) != null && StsKeywordIO.group.equals(group);
        }
    }

    class OverrideFilterDialog extends JDialog
    {
        StsJPanel panel = StsJPanel.addInsets();
        GridBagLayout gridBagLayout1 = new GridBagLayout();
        UserFilenameFilter filter = null;
        Frame frame = null;

        StsGroupBox defineBox = new StsGroupBox("Define Strings in Filename");
        StsIntFieldBean numberFieldsBean = null;
        StsStringFieldBean groupBean = null;
        StsStringFieldBean formatsBean = null;
        StsStringFieldBean tdBean = null;
        StsStringFieldBean markerBean = null;
        StsStringFieldBean logBean = null;
        StsStringFieldBean delimiterBean = null;

        String groupString = "well-dev";
        String formatsString = "las";
        String mkrString = "ref";
        String perfString = "perf";
        String tdString = "td";
        String logString = "logs";
        String delimiterString = ".";

        StsGroupBox orderBox = new StsGroupBox("Location of Strings in Filename");
        StsIntFieldBean grpOrderBean = null;
        StsIntFieldBean fmtOrderBean = null;
        StsIntFieldBean nameOrderBean = null;
        int numberFields = 3;
        int grpOrder = StsKeywordIO.GROUP;
        int fmtOrder = StsKeywordIO.FORMAT;
        int nameOrder = StsKeywordIO.NAME;

        StsButton okBtn = new StsButton();
        StsButton cancelBtn = new StsButton();

        public OverrideFilterDialog(Frame frame)
        {
            super(frame, "File Filter Definition", true);
            this.setLocation(frame.getLocation());
            this.frame = frame;
            try
            {
                initializeFields();
                constructBeans();
                jbInit();
                pack();
            }
            catch(Exception ex)
            {
                ex.printStackTrace();
            }
        }

        private void initializeFields()
        {
            numberFields = StsKeywordIO.getNumberFields();
            delimiterString = StsKeywordIO.getDelimiter();
            grpOrder = StsKeywordIO.getTokenOrder(StsKeywordIO.GROUP);
            fmtOrder = StsKeywordIO.getTokenOrder(StsKeywordIO.FORMAT);
            nameOrder = StsKeywordIO.getTokenOrder(StsKeywordIO.NAME);
            groupString = StsWellImport.getFilePrefix(StsWellImport.DEV);
            mkrString = StsWellImport.getFilePrefix(StsWellImport.REF);
            perfString = StsWellImport.getFilePrefix(StsWellImport.PERF);
            tdString = StsWellImport.getFilePrefix(StsWellImport.TD);
            logString = StsWellImport.getFilePrefix(StsWellImport.LOG);
        }

        private void constructBeans()
        {
            numberFieldsBean = new StsIntFieldBean(this, "numberFields", true, "Number of Fields:", true);
            groupBean = new StsStringFieldBean(this, "groupString", true, "Well Path:");
            formatsBean = new StsStringFieldBean(this, "formatsString", true, "Formats:");
            markerBean = new StsStringFieldBean(this, "mkrString", true, "Markers:");
            tdBean = new StsStringFieldBean(this, "tdString", true, "Time-Depth:");
            logBean = new StsStringFieldBean(this, "logString", true, "Logs:");
            delimiterBean = new StsStringFieldBean(this, "delimiterString", true, "Delimiter:");

            grpOrderBean = new StsIntFieldBean(this, "grpOrder", 1, 7, "Path/Log/Mkr:", true);
            fmtOrderBean = new StsIntFieldBean(this, "fmtOrder", 1, 7, "Format:", true);
            nameOrderBean = new StsIntFieldBean(this, "nameOrder", 1, 7, "Name:", true);

            okBtn = new StsButton("Ok", "Create well using input values", this, "okProcess", null);
            cancelBtn = new StsButton("Cancel", "Cancel well creation", this, "cancelProcess", null);
        }

        private void jbInit() throws Exception
        {
            panel.gbc.fill = panel.gbc.HORIZONTAL;
            panel.setBorder(BorderFactory.createEtchedBorder());
            panel.setLayout(gridBagLayout1);

            panel.gbc.gridwidth = 2;
            defineBox.addEndRow(numberFieldsBean);
            defineBox.addEndRow(groupBean);
            defineBox.addEndRow(formatsBean);
            defineBox.addEndRow(markerBean);
            defineBox.addEndRow(tdBean);
            defineBox.addEndRow(logBean);
            defineBox.addEndRow(delimiterBean);
            panel.add(defineBox);

            orderBox.addEndRow(grpOrderBean);
            orderBox.addEndRow(fmtOrderBean);
            orderBox.addEndRow(nameOrderBean);
            panel.add(orderBox);

            panel.gbc.gridwidth = 1;
            panel.addToRow(okBtn);
            panel.addEndRow(cancelBtn);

            this.getContentPane().add(panel);
        }
        public int getNumberFields() { return numberFields; }
        public String getGroupString() { return groupString; }
        public String getFormatsString() { return formatsString; }
        public int getGrpOrder() { return grpOrder + 1; }
        public int getFmtOrder() { return fmtOrder + 1; }
        public int getNameOrder() { return nameOrder + 1; }
        public String getTdString() { return tdString; }
        public String getLogString() { return logString; }
        public String getMkrString() { return mkrString; }
        public String getDelimiterString() { return delimiterString; }

        public void setNumberFields(int value) { numberFields = value; }
        public void setGroupString(String value) { groupString = value; }
        public void setFormatsString(String value) { formatsString = value; }
        public void setFmtOrder(int value) { fmtOrder = value - 1; }
        public void setGrpOrder(int value) { grpOrder = value - 1; }
        public void setNameOrder(int value) { nameOrder = value - 1; }
        public void setTdString(String value) { tdString = value; }
        public void setLogString(String value) { logString = value; }
        public void setMkrString(String value) { mkrString = value; }
        public void setDelimiterString(String value) { delimiterString = value; }

        public void okProcess()
        {
            // Parse the formats string
            if(formatsString == null)
                filter  = new UserFilenameFilter(getGroupString(), formats);
            else
            {
                int[] order = new int[] {grpOrder, fmtOrder, nameOrder};
                StsKeywordIO.setTokenOrder(order);
                String[] userFormats = new String[10];
                StringTokenizer stok = new StringTokenizer(formatsString,",");
                int nTokens = stok.countTokens();
                userFormats = new String[nTokens];
                for(int i=0; i<nTokens; i++)
                    userFormats[i] = stok.nextToken().trim();
                filter = new UserFilenameFilter(getGroupString(), userFormats);
                StsWellImport.setFilePrefixs(new String[] {groupString,logString,tdString,mkrString});
                StsKeywordIO.setDelimiter(delimiterString);
                StsKeywordIO.setNumberFields(numberFields);
            }
            this.hide();
        }

        public void cancelProcess()
        {
            this.setVisible(false);
        }

        public UserFilenameFilter getFilter() { return filter; }
    }

    class UserFilenameFilter extends StsFilenameFilter
    {
        String group = null;
        String[] formats = {"las"};
        public UserFilenameFilter(String group, String[] formats)  { this.group = group; this.formats = formats; }

        public boolean accept(File dir, String name)
        {
            if(StsWellImport.getCurrentDirectory() == null)
            {
                StsWellImport.setCurrentDirectory(dir.getAbsolutePath());
                wizard.constructAvailableFiles();
            }
            else
            {
                if(!StsWellImport.getCurrentDirectory().equals(dir.getAbsolutePath()))
                {
                    StsWellImport.setCurrentDirectory(dir.getAbsolutePath());
                    wizard.constructAvailableFiles();
                }
            }
            for(int i=0; i<formats.length; i++)
            {
                if ((name.endsWith("." + formats[i])) && (name.indexOf(group) != -1))
                    return true;
            }
            return false;
        }

        public void parseFilename(String filename) {}

        public String getFilenameName(String filename)
        {
            StsKeywordIO.parseAsciiFilename(filename);
            return StsKeywordIO.name;
        }
    }
}
