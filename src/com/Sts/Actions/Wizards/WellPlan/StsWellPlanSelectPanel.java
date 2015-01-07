package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsWellPlanSelectPanel extends JPanel implements ActionListener
{
    private JLabel jLabel1 = new JLabel();
    private JPanel selectedWellInfoPanel = new JPanel();
    private JPanel selectionPanel = new JPanel();
    private JPanel transferPanel = new JPanel();
    private JButton addBtn = new JButton();
    private JButton removeBtn = new JButton();
    private JButton viewWellBtn = new JButton();

    StsComboBoxFieldBean hUnitsBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean vUnitsBean = new StsComboBoxFieldBean();

    private JTextField currentDirectoryField = new JTextField();
    private JButton directoryBrowseButton = new JButton();

    private DefaultListModel selectedWellsListModel = new DefaultListModel();
    private DefaultListModel availableWellsListModel = new DefaultListModel();
    private String[] selectedWellPlanNames = null;

    static private JFileChooser chooseDirectory = null;

    private JScrollPane availableWellsScrollPane = new JScrollPane();
    private JList availableWellsList = new JList();
    private JScrollPane selectedWellsScrollPane = new JScrollPane();
    private JList selectedWellsList = new JList();

    private StsLoadWellPlanWizard wizard;
    private StsWellPlanSelect wizardStep;

    private byte vUnits, hUnits;
    private float vScalar = 1.0f, hScalar = 1.0f;

    public static int S2S = 0;
    public static int LAS = 1;
    public static int UT = 2;
   
    JRadioButton jRadioButton2 = new JRadioButton();
    ButtonGroup formatGroup = new ButtonGroup();
    JPanel filePanel = new JPanel();
//    VerticalFlowLayout verticalFlowLayout1 = new VerticalFlowLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    JCheckBox originalChk = new JCheckBox();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JPanel jPanel1 = new JPanel();
    GridBagLayout gridBagLayout4 = new GridBagLayout();
//    JTextField vScalarTxt = new JTextField();
//    JTextField hScalarTxt = new JTextField();
    GridBagLayout gridBagLayout5 = new GridBagLayout();

    public StsWellPlanSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsLoadWellPlanWizard)wizard;
        this.wizardStep = (StsWellPlanSelect)wizardStep;

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
        currentDirectoryField.setText(wizard.getModel().getProject().getRootDirString());
        availableWellsList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        hUnits = wizard.hUnits;
        vUnits = wizard.vUnits;
        hUnitsBean.initialize(this, "horzUnitsString", "Horizontal Units:  ", StsParameters.DIST_STRINGS);
        hUnitsBean.setSelectedItem(StsParameters.DIST_STRINGS[hUnits]);
        vUnitsBean.initialize(this, "vertUnitsString", "Vertical Units:  ", StsParameters.DIST_STRINGS);
        vUnitsBean.setSelectedItem(StsParameters.DIST_STRINGS[vUnits]);

        setAvailableFiles();
    }

    private void setAvailableFiles()
     {
         // Clear the List
         availableWellsListModel.removeAllElements();
		 StsFilenameFilter filenameFilter = new StsFilenameFilter("wellPlan", "obj");
		 String currentDirectory = currentDirectoryField.getText();
		 StsFileSet fileSet = StsFileSet.constructor(currentDirectory, filenameFilter);
		 StsAbstractFile[] files = fileSet.getFiles();
		 for(int n = 0; n < files.length; n++)
			 availableWellsListModel.addElement(files[n]);
   }
//
//  Get a list of all selected files...
//
//    public StsFile[] getSelectedFiles() { return selectedFiles; }

    public StsFile[] getSelectedFiles()
    {
      StsFile[] selectedFiles = new StsFile[selectedWellsListModel.getSize()];
      for (int i=0; i<selectedWellsListModel.getSize(); i++)
      {
           selectedFiles[i] = (StsFile)selectedWellsListModel.getElementAt(i);
      }
      return selectedFiles;
    }

// Get a list of all selected well objects...
//
/*
    public Object[] getSelectedWellObjects()
    {
      Object[] selectedWellObjects = null;

      selectedWellObjects = new Object[selectedWellsListModel.getSize()];
      for (int i=0; i<selectedWellsListModel.getSize(); i++)
        selectedWellObjects[i] = selectedWellsListModel.getElementAt(i);

      return selectedWellObjects;
    }
*/
    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout4);
        this.setMinimumSize(new Dimension(600, 600));
        jLabel1.setHorizontalAlignment(SwingConstants.CENTER);
        jLabel1.setText("Selected Well Information");
        selectedWellInfoPanel.setBorder(BorderFactory.createEtchedBorder());
        selectedWellInfoPanel.setLayout(gridBagLayout3);
        selectionPanel.setBorder(BorderFactory.createEtchedBorder());
        selectionPanel.setLayout(gridBagLayout2);
        this.setMinimumSize(new Dimension(0, 0));
        this.setPreferredSize(new Dimension(0, 0));
        removeBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        removeBtn.setMaximumSize(new Dimension(100, 20));
        removeBtn.setMinimumSize(new Dimension(100, 20));
        removeBtn.setPreferredSize(new Dimension(100, 20));
        removeBtn.setMargin(new Insets(0, 0, 0, 0));
        removeBtn.setText("< Remove");
        removeBtn.addActionListener(this);
        viewWellBtn.addActionListener(this);
        viewWellBtn.setText("View File");
        viewWellBtn.setMargin(new Insets(0, 0, 0, 0));
        viewWellBtn.setPreferredSize(new Dimension(100, 20));
        viewWellBtn.setMinimumSize(new Dimension(100, 20));
        viewWellBtn.setMaximumSize(new Dimension(100, 20));
        viewWellBtn.setFont(new java.awt.Font("Dialog", 0, 11));
//        transferPanel.setLayout(verticalFlowLayout1);
        transferPanel.setMaximumSize(new Dimension(50, 200));
        transferPanel.setMinimumSize(new Dimension(25, 140));
        transferPanel.setPreferredSize(new Dimension(50, 150));
        currentDirectoryField.setPreferredSize(new Dimension(100, 20));
//        currentDirectoryField.setText("jTextField1");
        ImageIcon icon = StsIcon.createIcon("dir16x32.gif");
        directoryBrowseButton.setIcon(icon);
        directoryBrowseButton.setMargin(new Insets(0, 0, 0, 0));
        directoryBrowseButton.addActionListener(this);
//        selectedWellsList.setMinimumSize(new Dimension(75, 200));
        selectedWellsList.setModel(selectedWellsListModel);
        availableWellsList.setModel(availableWellsListModel);
        addBtn.setText("Add >");
        addBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addBtn.setMaximumSize(new Dimension(100, 20));
        addBtn.setMinimumSize(new Dimension(100, 20));
        addBtn.setPreferredSize(new Dimension(100, 20));
        addBtn.setMargin(new Insets(0, 0, 0, 0));
        addBtn.addActionListener(this);
//        availableWellsScrollPane.setPreferredSize(new Dimension(150, 110));
//        selectedWellsScrollPane.setPreferredSize(new Dimension(100, 150));
        filePanel.setBorder(BorderFactory.createEtchedBorder());
        filePanel.setLayout(gridBagLayout1);
        originalChk.setText("Load from Original ASCII Files");
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout5);
        transferPanel.add(addBtn, null);
        transferPanel.add(removeBtn, null);
        transferPanel.add(viewWellBtn, null);

        originalChk.addActionListener(this);

        selectionPanel.add(selectedWellsScrollPane,  new GridBagConstraints(2, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 5, 5), 0, 0));
        selectionPanel.add(transferPanel,  new GridBagConstraints(1, 1, 1, 1, 0.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(0, 0, 5, 0), 100, 0));
        selectedWellsScrollPane.getViewport().add(selectedWellsList, null);
        selectionPanel.add(availableWellsScrollPane,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 5, 5), 0, 0));
        availableWellsScrollPane.getViewport().add(availableWellsList, null);
        selectionPanel.add(filePanel,  new GridBagConstraints(0, 0, 3, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 0, 0, 5), 0, -8));
        filePanel.add(directoryBrowseButton,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(7, 2, 0, 0), 0, 0));
        filePanel.add(currentDirectoryField,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(7, 0, 0, 5), 217, 0));
        filePanel.add(originalChk,   new GridBagConstraints(0, 1, 2, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 2, 5, 5), 0, 0));
        this.add(jPanel1,   new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 4, 0, 4), 373, 0));

        jPanel1.add(vUnitsBean,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 10, 0, 0), 0, 0));
        jPanel1.add(hUnitsBean,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 0, 0));

        this.add(selectedWellInfoPanel,  new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 4, 11, 4), 1, -50));
        selectedWellInfoPanel.add(jLabel1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 14, 91, 19), 219, 0));
        this.add(selectionPanel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 4, 0, 4), 2, 0));
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        int[] selectedIndices = null;
        String item = null;

        Object source = e.getSource();

        // select new directory
        if(source == directoryBrowseButton)
        {
            if (chooseDirectory == null) initializeChooseDirectory();
            chooseDirectory.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            chooseDirectory.setDialogTitle("Select or Enter Desired Directory and Press Open");
            chooseDirectory.setApproveButtonText("Open Directory");
            while(true)
            {
                chooseDirectory.showOpenDialog(null);
                File newDirectory = chooseDirectory.getSelectedFile();
                if(newDirectory == null)
                {
                    new StsMessage(wizard.getModel().win3d,  StsMessage.WARNING,  "No directory selected.");
                    break;
                }
                if (newDirectory.isDirectory())
                {
                    setCurrentDirectory(newDirectory.getAbsolutePath());
                    break;
                }
                else
                {
                    // File or nothing selected, strip off file name
                    String dirString = newDirectory.getPath();
                    newDirectory = new File(dirString.substring(0,dirString.lastIndexOf(File.separator)));
                    if(newDirectory.isDirectory())
                    {
                        setCurrentDirectory(newDirectory.getAbsolutePath());
                        break;
                    }
                    if(!StsYesNoDialog.questionValue(this,"Must select the directory that\n contains the Well Files.\n\n Continue?"))
                        break;
                }
            }
        }
        // Add selected volumes
        else if(source == addBtn)
        {
            selectedIndices = availableWellsList.getSelectedIndices();
            for(i=0; i< selectedIndices.length; i++)
            {
                if(selectedWellsListModel.indexOf(availableWellsList.getModel().getElementAt(selectedIndices[i])) == -1)
                   selectedWellsListModel.addElement(availableWellsList.getModel().getElementAt(selectedIndices[i]));
            }
        }
        // Remove selected volume(s)
        else if (source == removeBtn)
        {
            selectedIndices = selectedWellsList.getSelectedIndices();
            for(i= selectedIndices.length - 1; i>=0; i--)
                selectedWellsListModel.removeElementAt(selectedIndices[i]);
        }
        // Reset Available Files
        else if (source == originalChk)
            setAvailableFiles();

        else if(source == viewWellBtn)
        {
           String line;
           if(getCurrentFile().getName() == null)
           {
               new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "No Well Selected in Available List");
               return;
           }

           StsFileViewDialog dialog = new StsFileViewDialog(wizard.frame,"Well File View", false);
           dialog.setVisible(true);

           dialog.setViewTitle("Well File - " + getCurrentFile().getName());
           try
           {
               BufferedReader wellFile = new BufferedReader(new FileReader(getCurrentFile()));
               while (true) {
                   line = wellFile.readLine();
                   if (line == null)
                       break;
                   dialog.appendLine(line);
               }
           }
           catch (Exception ex) {
               StsMessageFiles.infoMessage("Unable to view well file.");
           }
       }
    }
/*
    public Object[] getSelectedWellObjects()
    {
       Object[] selectedWellObjects = null;

       selectedWellObjects = new Object[selectedWellsListModel.getSize()];
       for (int i=0; i<selectedWellsListModel.getSize(); i++)
         selectedWellObjects[i] = selectedWellsListModel.getElementAt(i);

       return selectedWellObjects;
     }
*/
    private boolean isInSelectedList(String item)
    {
        return true;
    }

    private void initializeChooseDirectory()
    {
        chooseDirectory = new JFileChooser(StsWellImport.getCurrentDirectory());
    }

    private void setCurrentDirectory(String directory)
    {
        StsWellImport.setCurrentDirectory(directory);
        currentDirectoryField.setText(directory);
        setAvailableFiles();
    }

    public boolean getReloadASCII()
    {
        return originalChk.isSelected();
    }
    public void setReloadASCII(boolean value)
    {
        originalChk.setSelected(value);
    }

    public File getCurrentFile()
    {
        String filename = StsWellImport.getCurrentDirectory() + File.separator + StsWellImport.getWellFilename(availableWellsList.getSelectedIndex());
        return new File(filename);
    }

    public String getHorzUnitsString() { return StsParameters.DIST_STRINGS[hUnits]; }
    public String getVertUnitsString() { return StsParameters.DIST_STRINGS[vUnits]; }
    public byte getHorzUnits() { return hUnits; }
    public byte getVertUnits() { return vUnits; }
    public void setHorzUnitsString(String unitString)
    {
        hUnits = StsParameters.getDistanceUnitsFromString(unitString);
//        hScalar = wizard.getModel().getProject().getXyScalar(hUnits);
//        hScalarTxt.setText(Float.toString(hScalar));
    }
//    public float getHorzScalar() { return hScalar; }

    public void setVertUnitsString(String unitString)
    {
        vUnits = StsParameters.getDistanceUnitsFromString(unitString);
//        vScalar = wizard.getModel().getProject().getDepthScalar(vUnits);
//        vScalarTxt.setText(Float.toString(vScalar));
    }
//    public float getVertScalar() { return vScalar; }
}


