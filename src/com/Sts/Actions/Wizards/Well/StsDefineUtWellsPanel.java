package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.event.*;
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

public class StsDefineUtWellsPanel extends JPanel implements ListSelectionListener, ActionListener
{
    private StsWizard wizard;
    private StsDefineUtWells wizardStep;

    private StsPoint[] points = null;
    private int currentIdx = -1;
    private byte[] type = null;
    private int xOriginCol = 9;
    private int yOriginCol = 10;
    private int zOriginCol = -1;
    private double zValue = 0.0f;
    private float scalar = 1.0f;

    StsGroupBox propertiesBox = new StsGroupBox("Top Hole Location");
    StsDoubleFieldBean xBean = new StsDoubleFieldBean();
    StsDoubleFieldBean yBean = new StsDoubleFieldBean();
    StsDoubleFieldBean zBean = new StsDoubleFieldBean();
    JButton defineBtn = new JButton("Define From File");
    StsComboBoxFieldBean typeListBean = new StsComboBoxFieldBean(this, "type", "Algorithm: ");

    JPanel jPanel1 = new JPanel();
    JScrollPane jScrollPane1 = new JScrollPane();
    JList fileList = new JList();
    JLabel jLabel1 = new JLabel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JDialog dialog = null;
    StsWellDefinePanel definePanel = null;
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    public StsDefineUtWellsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = (StsDefineUtWells) wizardStep;

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
        String[] selectedWellnames = wizardStep.getSelectedWells();
        String[] selectedWellFilenames = StsWellImport.getSelectedWellFilenames(selectedWellnames, StsWellImport.UTFILES, true);
        int nWells = selectedWellFilenames.length;
        points = new StsPoint[nWells];
        type = new byte[nWells];

        double[] origin =
            { (double) wizard.getModel().getProject().getXOrigin(),
            (double) wizard.getModel().getProject().getYOrigin(), 0.0f};

        for (int i = 0; i < nWells; i++)
        {
            points[i] = new StsPoint(origin);
            type[i] = StsUTKeywordIO.MINCURVE;
        }

        fileList.setListData(selectedWellFilenames);
        typeListBean.setListItems(StsUTKeywordIO.algorithms);
    }

    void jbInit() throws Exception
    {
        xBean.initialize(this, "Vx", true, "Slot X:");
        yBean.initialize(this, "Vy", true, "Slot Y:");
        zBean.initialize(this, "Vz", true, "Kelly Bushing Elevation:");
        defineBtn.addActionListener(this);

        propertiesBox.add(xBean);
        propertiesBox.add(yBean);
        propertiesBox.add(zBean);
        propertiesBox.add(typeListBean);

        fileList.addListSelectionListener(this);

        jPanel1.setLayout(gridBagLayout3);
        this.setLayout(gridBagLayout2);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jLabel1.setText("Selected UT Files");
        this.add(jPanel1, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
                                                 , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 3, 2, 4), 0, 3));
        jPanel1.add(jScrollPane1, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(9, 0, 0, 7), 51, 48));
        jPanel1.add(jLabel1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 0, 0, 202), 33, 0));
        jPanel1.add(propertiesBox, new GridBagConstraints(0, 3, 1, 1, 1.0, 1.0
            , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 2, 0), 182, 48));
        propertiesBox.add(defineBtn);
        jScrollPane1.getViewport().add(fileList, null);
    }

    public void valueChanged(ListSelectionEvent e)
    {
        Object source = e.getSource();

        if (source == fileList)
        {
            if (fileList.getSelectedIndex() != -1)
            {
                currentIdx = fileList.getSelectedIndex();
                showTopHoleLocation();
            }
        }
    }

    public String getSelectedFile()
    {
        return (String) fileList.getSelectedValue();
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        int[] selectedIndices = null;

        Object source = e.getSource();

        if (source == defineBtn)
        {
            if (getSelectedFile() == null)
            {
                new StsMessage(wizard.getModel().win3d, StsMessage.ERROR, "Must select a well file first.");
                return;
            }
            // Show Trace Display Dialog
            dialog = new JDialog(wizard.frame, true);
            definePanel = new StsWellDefinePanel(wizard, wizardStep);
            definePanel.initialize();
            dialog.setSize(700, 400);
            dialog.setTitle("Well File Define / View");
            dialog.getContentPane().add(definePanel);
            dialog.setVisible(true);
        }
    }

    public void exitFileView()
    {
        xOriginCol = definePanel.getXCol();
        yOriginCol = definePanel.getYCol();
        zValue = definePanel.getZ();
        if ( (xOriginCol != -1) && (yOriginCol != -1))
        {
            int[] indices = fileList.getSelectedIndices();
            for (int i = 0; i < indices.length; i++)
            {
                fileList.setSelectedIndex(indices[i]);
                setVz(zValue);
                if (!getTopLocationsFromFile( (String) fileList.getModel().getElementAt(indices[i])))
                {
                    new StsMessage(wizard.getModel().win3d, StsMessage.WARNING,
                                   "Unable to extract origin from file: " + (String) fileList.getModel().getElementAt(indices[i]));
                }
                else
                {
                    showTopHoleLocation();
                }
            }
        }
        dialog.setVisible(false);
        dialog.dispose();
    }

    public boolean getTopLocationsFromFile(String filename)
    {
        String[] tokens = null;
        double xValue = 0.0f, yValue = 0.0f;

        if (filename == null)
        {
            return false;
        }

        String dir = StsWellImport.getCurrentDirectory() + File.separator;
        StsFile file = StsFile.constructor(dir, filename);
        StsAsciiFile asciiFile = new StsAsciiFile(file);
        if(!asciiFile.openReadWithErrorMessage()) return false;

        try
        {
            while (true)
            {
                tokens = asciiFile.getTokens();
                try
                {
                    xValue = Double.parseDouble(tokens[0]);
                    break;
                }
                catch (Exception e)
                {
                    continue;
                }
            }
            int nLines = 0;
            int numCols = tokens.length;
            if ( (numCols < xOriginCol) || (numCols < yOriginCol) || (numCols < zOriginCol))
            {
                return false;
            }

            while (tokens != null)
            {
                if (++nLines > 100)
                {
                    break;
                }
                int nTokens = tokens.length;
                try
                {
                    xValue = Double.parseDouble(tokens[xOriginCol]);
                    yValue = Double.parseDouble(tokens[yOriginCol]);
                    setVx(xValue);
                    setVy(yValue);
                    break;
                }
                catch (Exception e)
                {
                    tokens = asciiFile.getTokens();
                    continue;
                }
            }

        }
        catch (Exception e)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "File read error for " +
                           filename + ": " + e.getMessage());
            return false;
        }

        return true;
    }

    private void showTopHoleLocation()
    {
        xBean.setValue(points[currentIdx].getX());
        yBean.setValue(points[currentIdx].getY());
        zBean.setValue(points[currentIdx].getZ());
    }

    public void setVx(double x)
    {
        int[] selectedIndices = fileList.getSelectedIndices();
        for (int i = 0; i < fileList.getSelectedIndices().length; i++)
        {
            points[selectedIndices[i]].setX( (float) x);
        }
    }

    public void setVy(double y)
    {
        int[] selectedIndices = fileList.getSelectedIndices();
        for (int i = 0; i < fileList.getSelectedIndices().length; i++)
        {
            points[selectedIndices[i]].setY( (float) y);
        }
    }

    public double getVx()
    {
        if (currentIdx != -1)
        {
            return points[currentIdx].getX();
        }
        else
        {
            return -1.0f;
        }
    }

    public double getVy()
    {
        if (currentIdx != -1)
        {
            return points[currentIdx].getY();
        }
        else
        {
            return -1.0f;
        }
    }

    public void setVz(double z)
    {
        int[] selectedIndices = fileList.getSelectedIndices();
        for (int i = 0; i < fileList.getSelectedIndices().length; i++)
        {
            if (z < 0.0) // Already in depth
            {
                points[selectedIndices[i]].setZ( (float) z);
            }
            else // Convert elevation to depth
            {
                points[selectedIndices[i]].setZ( (float) - z);
            }
        }
    }

    public double getVz()
    {
        if (currentIdx != -1)
        {
            return points[currentIdx].getZ();
        }
        else
        {
            return -1.0f;
        }
    }

    public StsPoint[] getTopHoleLocations()
    {
        return points;
    }

    public void setType(String type)
    {
        int[] selectedIndices = fileList.getSelectedIndices();
        for (int i = 0; i < fileList.getSelectedIndices().length; i++)
        {
            this.type[selectedIndices[i]] = StsUTKeywordIO.getAlgorithmTypeFromName(type);
        }
    }

    public String getType()
    {
        if (currentIdx >= 0)
        {
            return StsUTKeywordIO.getAlgorithmNameFromType(type[currentIdx]);
        }
        else
        {
            return "Unknown";
        }
    }

    public byte[] getTypes()
    {
        return type;
    }

}
