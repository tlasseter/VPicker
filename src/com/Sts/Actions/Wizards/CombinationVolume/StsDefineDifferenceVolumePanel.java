package com.Sts.Actions.Wizards.CombinationVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Histogram.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineDifferenceVolumePanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;

    private StsSeismicVolume[] seismicVolumes = null;
    private StsVirtualVolume[] virtualVolumes = null;

    private int nSeismicVolumes = 0;
    private int nVirtualVolumes = 0;

    private int selectedSeismicOneIndex = 0;
    private int selectedSeismicTwoIndex = 0;

    private StsFloatFieldBean dataMinBean;
    private StsFloatFieldBean dataMaxBean;
    private float dataMin = StsParameters.largeFloat;
    private float dataMax = StsParameters.smallFloat;

    JComboBox seismicOneCombo = new JComboBox();
    JLabel seismicOneLabel = new JLabel();
    StsSeismicVolume volOne = null;
    JComboBox seismicTwoCombo = new JComboBox();
    JLabel seismicTwoLabel = new JLabel();
    StsSeismicVolume volTwo = null;

    JPanel jPanel3 = new JPanel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JLabel seismicTwoLabel1 = new JLabel();
    JLabel seismicTwoLabel2 = new JLabel();
    JLabel maxLabel = new JLabel();
    JPanel jPanel1 = new JPanel();
    JLabel minLabel = new JLabel();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    StsHistogramPanel histogram = new StsHistogramPanel(StsHistogramPanel.HORIZONTAL);
    JLabel histogramMaxLabel = new JLabel("Maximum");
    JLabel histogramMinLabel = new JLabel("Minimum");
    JLabel histogramMaxPctLabel = new JLabel("Clip %");
    JLabel histogramMinPctLabel = new JLabel("Clip %");
    JPanel dataHistogram = new JPanel();

    private DecimalFormat numberFormat = new DecimalFormat("##0.0#");

    public StsDefineDifferenceVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
            dataMinBean = new StsFloatFieldBean(this, "dataMin");
            dataMaxBean = new StsFloatFieldBean(this, "dataMax");
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
        //
        // Add all existing seismic volumes
        //
        seismicVolumes = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
        nSeismicVolumes = seismicVolumes.length;
        for(int v = 0; v < nSeismicVolumes; v++)
        {
            seismicOneCombo.addItem(seismicVolumes[v].getName());
            seismicTwoCombo.addItem(seismicVolumes[v].getName());
        }
        virtualVolumes = ((StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class)).getVirtualVolumes();
        nVirtualVolumes = virtualVolumes.length;
        for(int v = 0; v < nVirtualVolumes; v++)
        {
            seismicOneCombo.addItem(virtualVolumes[v].getName());
            seismicTwoCombo.addItem(virtualVolumes[v].getName());
        }
        if((nSeismicVolumes + nVirtualVolumes) > 0)
        {
            seismicOneCombo.setSelectedIndex(0);
            if((nSeismicVolumes + nVirtualVolumes) > 1)
                seismicTwoCombo.setSelectedIndex(1);
            else
                seismicTwoCombo.setSelectedIndex(0);
        }
        dataMax = Math.max(volOne.getDataMax(), volTwo.getDataMax());
        dataMin = Math.min(volOne.getDataMin(), volTwo.getDataMin());
        resetLimits();
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);
        seismicOneLabel.setFont(new java.awt.Font("Dialog", 3, 11));
        seismicOneLabel.setHorizontalAlignment(SwingConstants.RIGHT);
        seismicOneLabel.setText("PostStack3d One:");
        seismicTwoLabel.setFont(new java.awt.Font("Dialog", 3, 11));
        seismicTwoLabel.setHorizontalAlignment(SwingConstants.RIGHT);
        seismicTwoLabel.setText("PostStack3d Two:");

        seismicOneCombo.setToolTipText("Select first seismic volume");
        seismicTwoCombo.setToolTipText("Select optional second seismic volume or specify scalar");
        jPanel3.setBorder(BorderFactory.createEtchedBorder());
        jPanel3.setLayout(gridBagLayout3);
        seismicTwoLabel1.setText("minus");
        seismicTwoLabel1.setHorizontalAlignment(SwingConstants.CENTER);
        seismicTwoLabel1.setFont(new java.awt.Font("Dialog", 3, 11));
        seismicTwoLabel2.setText("Null Limit Based on PostStack3d One");
        seismicTwoLabel2.setHorizontalAlignment(SwingConstants.CENTER);
        seismicTwoLabel2.setFont(new java.awt.Font("Dialog", 3, 11));
        maxLabel.setFont(new java.awt.Font("Dialog", 0, 11));
        maxLabel.setText("Maximum:");
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout1);
        minLabel.setFont(new java.awt.Font("Dialog", 0, 11));
        minLabel.setText("Minimum:");

        histogramMinLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMinLabel.setForeground(Color.gray);
        histogramMinPctLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMinPctLabel.setForeground(Color.gray);
        histogramMaxLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMaxLabel.setForeground(Color.gray);
        histogramMaxPctLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMaxPctLabel.setForeground(Color.gray);

        jPanel3.add(seismicOneCombo,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(6, 0, 0, 3), 245, 4));
        jPanel3.add(seismicTwoCombo,  new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 3), 245, 4));
        jPanel3.add(seismicTwoLabel1,  new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 84, 0, 91), 68, 4));
        jPanel3.add(jPanel1,  new GridBagConstraints(0, 3, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(6, 0, 4, 3), 0, 0));
        this.add(jPanel3, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 3, 5, 4), 0, 0));

        jPanel1.add(seismicTwoLabel2, new GridBagConstraints(0, 0, 4, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(minLabel,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 2, 2, 2), 0, 0));
        jPanel1.add(dataMinBean,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(15, 2, 10, 0), 0, 0));
        jPanel1.add(maxLabel,  new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 3, 0, 0), 0, 0));
        jPanel1.add(dataMaxBean,  new GridBagConstraints(3, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(10, 3, 0, 0), 0, 0));
        jPanel1.add(dataHistogram,  new GridBagConstraints(0, 3, 4, 2, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 0, 0, 0), 0, 0));

        jPanel3.add(seismicTwoLabel,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 11, 0, 0), 33, 4));
        jPanel3.add(seismicOneLabel,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 3, 0, 0), 42, 3));

        seismicOneCombo.addActionListener(this);
        seismicTwoCombo.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        if(source == seismicOneCombo)
        {
            selectedSeismicOneIndex = seismicOneCombo.getSelectedIndex();
            volOne = getSelectedSeismicOneVolume();
        }
        else if(source == seismicTwoCombo)
        {
            selectedSeismicTwoIndex = seismicTwoCombo.getSelectedIndex();
            volTwo = getSelectedSeismicTwoVolume();
        }
        resetLimits();
    }

    private void resetLimits()
    {
        float max = 0.0f, min = 0.0f;

        if((volOne == null) || (volTwo == null))
            return;

        max = Math.max(volOne.getDataMax(), volTwo.getDataMax());
        min = Math.min(volOne.getDataMin(), volTwo.getDataMin());

        if((dataMin < min) || (dataMin == StsParameters.largeFloat))
            dataMin = min;
        if((dataMax > max) || (dataMax == StsParameters.smallFloat))
            dataMax = max;

        clearHistogramPanel();
        displayHistogram(volOne);

        dataMinBean.setValueAndRange(dataMin, min, max);
        dataMaxBean.setValueAndRange(dataMax, min, max);

        setDataMin(dataMin);
        setDataMax(dataMax);
    }

    public void displayHistogram(StsSeismicVolume volume)
    {
        histogram.updateData(volume.dataHist, volOne.getDataMin(), volOne.getDataMax());

        String label = new String(numberFormat.format(histogram.getBottomPercentageClipped()) + "%");
        histogramMinPctLabel.setText(label);
        dataHistogram.add(histogramMinPctLabel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 30, 0));

        dataHistogram.add(histogram,  new GridBagConstraints(1, 0, 3, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 200, 0));

        label = new String(" " + numberFormat.format(histogram.getTopPercentageClipped()) + "%");
        histogramMaxPctLabel.setText(label);
        dataHistogram.add(histogramMaxPctLabel,  new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 30, 0));

        histogram.repaint();
        validate();
    }

    public void clearHistogramPanel()
    {
        dataHistogram.remove(histogram);
        dataHistogram.remove(histogramMaxLabel);
        dataHistogram.remove(histogramMinLabel);
        dataHistogram.remove(histogramMinPctLabel);
        dataHistogram.remove(histogramMaxPctLabel);
        histogram.clearAll();
    }

    public StsSeismicVolume getSelectedSeismicOneVolume()
    {
        if(selectedSeismicOneIndex > seismicVolumes.length - 1)
            return virtualVolumes[selectedSeismicOneIndex - seismicVolumes.length];
        else
            return seismicVolumes[selectedSeismicOneIndex];
    }

    public StsSeismicVolume getSelectedSeismicTwoVolume()
    {
        if(selectedSeismicTwoIndex > seismicVolumes.length - 1)
            return virtualVolumes[selectedSeismicTwoIndex - seismicVolumes.length];
        else
            return seismicVolumes[selectedSeismicTwoIndex];
    }

    public float getDataMax() { return dataMax; }
    public void setDataMax(float value)
    {
        dataMax = value;
        if(histogram != null)
        {
            int minIdx = (int)((dataMin - volOne.getDataMin()) * (histogram.getNumberIndices()/(volOne.getDataMax() - volOne.getDataMin())));
            int maxIdx = (int)((dataMax - volOne.getDataMin()) * (histogram.getNumberIndices()/(volOne.getDataMax() - volOne.getDataMin())));
            histogram.setClip(minIdx, maxIdx);
            histogram.repaint();
            histogramMaxPctLabel.setText(new String(numberFormat.format(histogram.getTopPercentageClipped()) + "%"));
        }

    }
    public float getDataMin() { return dataMin; }
    public void setDataMin(float value)
    {
        dataMin = value;
        if(histogram != null)
        {
            int minIdx = (int)((dataMin - volOne.getDataMin()) * (histogram.getNumberIndices()/(volOne.getDataMax() - volOne.getDataMin())));
            int maxIdx = (int)((dataMax - volOne.getDataMin()) * (histogram.getNumberIndices()/(volOne.getDataMax() - volOne.getDataMin())));
            histogram.setClip(minIdx, maxIdx);
            histogram.repaint();
            histogramMinPctLabel.setText(new String(numberFormat.format(histogram.getBottomPercentageClipped()) + "%"));
        }

    }

}
