package com.Sts.Actions.Wizards.CrossPlot;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

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

public class StsDefineCrossplotPanel extends JPanel implements ActionListener, ItemListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;
    private StsSeismicVolume[] seismicVolumes = null;
    private int nSeismicVolumes = 0;
    private int[] axisVolumeIndices;
    private int attributeVolumeIdx = 0;
    private int nAxes = 2;
    private File polygonFile = null;

    private StsTypeLibrary[] typeLibraries = null;
    DefaultListModel libraryListModel = new DefaultListModel();
    private int nTypeLibraries = 0;
    private int selectedLibrary = 0;

    JTextField nameTextField = new JTextField();
    JLabel nameLabel = new JLabel();
    JButton polygonViewButton = new JButton();

    JComboBox[] axisVolumes = new JComboBox[] { new JComboBox(), new JComboBox() };
    JComboBox attributeVolume = new JComboBox();
    JLabel[] axisVolLabels = new JLabel[] { new JLabel(), new JLabel() };
    JLabel attributeLabel = new JLabel();
    JLabel polygonFileText = new JLabel();
/*
    JComboBox xAxisVolume = new JComboBox();
    JLabel xAxisVolLabel = new JLabel();
    JComboBox yAxisVolume = new JComboBox();
    JLabel yAxisVolLabel = new JLabel();
*/
    JButton importPolygonButton = new JButton();
    JComboBox typesLibraryCombo = new JComboBox();
    JLabel jLabel1 = new JLabel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsDefineCrossplotPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
            jbInit();
//            classInitialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
        seismicVolumes = ((StsCrossplotWizard)wizard).getSeismicVolumes();
        nSeismicVolumes = seismicVolumes.length;
        axisVolumeIndices = new int[nSeismicVolumes];

        for(int a = 0; a < nAxes; a++)
        {
            for(int v = 0; v < nSeismicVolumes; v++)
                axisVolumes[a].addItem(seismicVolumes[v].getName());

            setAxisVolumeSelected(a, a);
        }
        attributeVolume.addItem("None");
        for(int v = 0; v < nSeismicVolumes; v++)
        {
            attributeVolume.addItem(seismicVolumes[v].getName());
        }
        setAttributeVolumeSelected(0);

        typeLibraries = (StsTypeLibrary[])model.getCastObjectList(StsTypeLibrary.class);
        nTypeLibraries = typeLibraries.length;
        if(nTypeLibraries > 0)
        {
            for (int n = 0; n < nTypeLibraries; n++)
                typesLibraryCombo.addItem(typeLibraries[n].getName());
            typesLibraryCombo.setSelectedItem(StsTypeLibrary.genericLibraryName);
        }
    }

    private void checkSetAttributeVolume(int index)
    {
        // if this index is already assigned to this axis, return
        if(attributeVolumeIdx == index)
            return;

        // assign volumeIndex to new axis
        attributeVolumeIdx = index;
    }
    private void checkSetAxisVolumeSelected(int nAxis, int volumeIndex)
    {
        // if this index is already assigned to this axis, return
        if(axisVolumeIndices[nAxis] == volumeIndex) return;

        // remove this volumeIndex from current axis
        for(int a = 0; a < nAxes; a++)
            if(axisVolumeIndices[a] == volumeIndex)
                 axisVolumeIndices[a] = -1;

        // assign volumeIndex to new axis
        axisVolumeIndices[nAxis] = volumeIndex;

        // get list of unassigned volumes
        int nUnassigned = 0;
        boolean [] unassigned = new boolean[nSeismicVolumes];
        for(int v = 0; v < nSeismicVolumes; v++)
            unassigned[v] = true;
        for(int a = 0; a < nAxes; a++)
        {
            volumeIndex = axisVolumeIndices[a];
            if(volumeIndex != -1) unassigned[volumeIndex] = false;
        }

        // for axes without volume assigned, assign one
        for(int a = 0; a < nAxes; a++)
        {
            if(axisVolumeIndices[a] == -1)
            {
                for(int v = 0; v < nSeismicVolumes; v++)
                    if(unassigned[v])
                    {
                        setAxisVolumeSelected(a, v);
                        unassigned[v] = false;
                    }
            }
        }
    }

    private void setAxisVolumeSelected(int nAxis, int volumeIndex)
    {
        axisVolumes[nAxis].setSelectedIndex(volumeIndex);
    }
    private void setAttributeVolumeSelected(int volumeIndex)
    {
        attributeVolume.setSelectedIndex(volumeIndex);
    }

    public int getAttributeVolumeIndex() { return attributeVolumeIdx; }
    public int[] getAxisVolumeIndices()
    {
        return axisVolumeIndices;
    }

    void jbInit() throws Exception
    {
        nameTextField.setText("");
        this.setLayout(gridBagLayout1);
        nameLabel.setText("Crossplot Name:");
        axisVolLabels[0].setText("X-axis PostStack3d:");
        axisVolLabels[1].setText("Y-axis PostStack3d:");
        attributeLabel.setText("Attribute PostStack3d:");
        importPolygonButton.setText("Import Polygons");
        importPolygonButton.setToolTipText("Import polygons from ASCII file to apply to this crossplot");
        polygonFileText.setText("none");
        polygonFileText.setBorder(BorderFactory.createEtchedBorder());
        polygonViewButton.setText("View");
        polygonViewButton.setToolTipText("View the selected polygon file");

        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setText("Types Library:");
        typesLibraryCombo.addItemListener(this);
        typesLibraryCombo.setToolTipText("Select the colors/lithologic types that will be available for polygon definition");

        this.add(nameTextField,  new GridBagConstraints(1, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 5), 0, 0));
        this.add(nameLabel,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 2), 0, 0));

        this.add(axisVolumes[0],  new GridBagConstraints(1, 1, 2, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 5), 0, 0));
        this.add(axisVolLabels[0],   new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 2), 0, 0));

        this.add(axisVolumes[1],  new GridBagConstraints(1, 2, 2, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 5), 0, 0));
        this.add(axisVolLabels[1],   new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 2), 0, 0));

        this.add(attributeVolume,  new GridBagConstraints(1, 3, 2, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 5), 0, 0));
        this.add(attributeLabel,   new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 2), 0, 0));
        if(((StsCrossplotWizard)wizard).getNumberOfAvailableVolumes() < 3)
        	attributeVolume.setEnabled(false);
        
        this.add(importPolygonButton,  new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 2), 0, 0));
        this.add(polygonFileText,  new GridBagConstraints(1, 4, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 2), 0, 0));
        this.add(polygonViewButton,  new GridBagConstraints(2, 4, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 5), 0, 0));

        this.add(typesLibraryCombo,  new GridBagConstraints(1, 5, 2, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 5), 0, 0));
        this.add(jLabel1,  new GridBagConstraints(0, 5, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 2), 0, 0));


        for(int n = 0; n < nAxes; n++)
            axisVolumes[n].addActionListener(this);

        attributeVolume.addActionListener(this);
        importPolygonButton.addActionListener(this);
        polygonViewButton.addActionListener(this);
    }

    public void itemStateChanged(ItemEvent e)
    {
        Object source = e.getSource();
        if(source == typesLibraryCombo)
        {
            ((StsCrossplotWizard)wizard).setSelectedLibrary(
                typeLibraries[typesLibraryCombo.getSelectedIndex()]);
        }
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        for(int nAxis = 0; nAxis < nAxes; nAxis++)
        {
            if(source == axisVolumes[nAxis])
            {
                int index = axisVolumes[nAxis].getSelectedIndex();
                checkSetAxisVolumeSelected(nAxis, index);
            }
        }
        if(source == attributeVolume)
        {
            int index = attributeVolume.getSelectedIndex();
            checkSetAttributeVolume(index);
        }
        else if(source == importPolygonButton)
        {
            JFileChooser chooser = new JFileChooser(model.getProject().getRootDirString());
            chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
            chooser.setDialogTitle("Select a polygon file");
            chooser.setApproveButtonText("Open File");
            while(true)
            {
                chooser.showOpenDialog(null);
                polygonFile = chooser.getSelectedFile();
                // Is the file a polygon file
				// jbw check for null
				if (polygonFile != null) {
					if(polygonFile.getPath().indexOf("polygons.txt") != -1)
					{
						polygonFileText.setText(polygonFile.getName());
						break;
					}
				}
            }
        }
        else if(source == polygonViewButton)
        {
           String line;

           if(getPolygonFile() == null)
           {
               new StsMessage(model.win3d, StsMessage.WARNING, "No Polygon File Specified.");
               return;
           }

           StsFileViewDialog dialog = new StsFileViewDialog(wizard.frame,"Polygon File View", false);
           dialog.setVisible(true);

           dialog.setViewTitle("Polygon File - " + polygonFile.getName());
           try
           {
               BufferedReader polygonFile = new BufferedReader(new FileReader(getPolygonFile()));
               while (true)
               {
                   line = polygonFile.readLine();
                   if (line == null)
                       break;
                   dialog.appendLine(line);
               }
           }
           catch (Exception ex)
           {
               StsMessageFiles.infoMessage("Unable to view polygon file.");
           }
       }
    }
    public File getPolygonFile() { return polygonFile; }

    public String getName() { return nameTextField.getText(); }

    public StsSeismicVolume[] getAxisSeismicVolumes()
    {
        StsSeismicVolume[] axisVolumes = new StsSeismicVolume[nAxes];
        for(int a = 0; a < nAxes; a++)
        {
            int selectedIndex = axisVolumeIndices[a];
            if(selectedIndex >= 0)
                axisVolumes[a] = seismicVolumes[selectedIndex];
        }
        return axisVolumes;
    }
    public StsSeismicVolume getAttributeSeismicVolume()
    {
        if(attributeVolumeIdx == 0) return null;
        else return seismicVolumes[attributeVolumeIdx-1];
}
}
