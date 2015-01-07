package com.Sts.Actions.Wizards.CombinationVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Histogram.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
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

public class StsDefineCummulativeDifferenceVolumePanel extends JPanel implements ActionListener, ListSelectionListener
{
    private StsWizard wizard;
    private StsDefineCummulativeDifferenceVolume wizardStep;

    private JPanel selectionPanel = new JPanel();
    private JPanel transferPanel = new JPanel();
    private JButton addBtn = new JButton();
    private JButton addAllBtn = new JButton();
    private JButton removeBtn = new JButton();
    private JButton removeAllBtn = new JButton();

    private StsSeismicVolume selectedVol = null;
    private int selectedIdx = 0;
    private DefaultListModel availableVolsListModel = new DefaultListModel();
    private JScrollPane availableVolsScrollPane = new JScrollPane();
    private JList availableVolsList = new JList();

    private DefaultListModel selectedVolsListModel = new DefaultListModel();
    private JScrollPane selectedVolsScrollPane = new JScrollPane();
    private JList selectedVolsList = new JList();

    static final int UP = 0;
    static final int DOWN = 1;

    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    private StsModel model = null;
    ButtonGroup modeGroup = new ButtonGroup();
    JButton upBtn = new JButton();
    JButton downBtn = new JButton();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JPanel histPanel = new JPanel();
    GridBagLayout gridBagLayout4 = new GridBagLayout();

    StsHistogramPanel histogram = new StsHistogramPanel(StsHistogramPanel.HORIZONTAL);
    JLabel histogramMaxLabel = new JLabel("Maximum");
    JLabel histogramMinLabel = new JLabel("Minimum");
    JLabel histogramMaxPctLabel = new JLabel("Clip %");
    JLabel histogramMinPctLabel = new JLabel("Clip %");
    JLabel maxLabel = new JLabel("Maximum");
    JLabel minLabel = new JLabel("Maximum");
    JLabel histLabel = new JLabel();
    JPanel dataHistogram = new JPanel();
    private StsFloatFieldBean dataMinBean;
    private StsFloatFieldBean dataMaxBean;
    private float dataMin = StsParameters.largeFloat;
    private float dataMax = StsParameters.smallFloat;
    public StsCheckbox autoChk;

    private DecimalFormat numberFormat = new DecimalFormat("##0.0#");

    public StsDefineCummulativeDifferenceVolumePanel(StsWizard wizard, StsDefineCummulativeDifferenceVolume wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
            dataMinBean = new StsFloatFieldBean(this, "dataMin");
            dataMaxBean = new StsFloatFieldBean(this, "dataMax");
            jbInit();
            initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
        StsSeismicVolume[] volumes = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);

        int nSeismicVols = volumes.length;
        availableVolsListModel.removeAllElements();
        for(int i=0; i<nSeismicVols; i++)
            availableVolsListModel.addElement(volumes[i].getName());

        volumes = ((StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class)).getVirtualVolumes();
        if(volumes != null)
        {
        int nVirtualVols = volumes.length;
        for(int i=0; i<nVirtualVols; i++)
            availableVolsListModel.addElement(volumes[i].getName());
        }
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout4);
        selectionPanel.setBorder(BorderFactory.createEtchedBorder());
        selectionPanel.setLayout(gridBagLayout2);
        removeAllBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        removeAllBtn.setMaximumSize(new Dimension(100, 20));
        removeAllBtn.setMinimumSize(new Dimension(100, 20));
        removeAllBtn.setPreferredSize(new Dimension(100, 20));
        removeAllBtn.setMargin(new Insets(0, 0, 0, 0));
        removeAllBtn.setText("<< Remove All");
        removeAllBtn.addActionListener(this);
        removeBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        removeBtn.setMaximumSize(new Dimension(100, 20));
        removeBtn.setMinimumSize(new Dimension(100, 20));
        removeBtn.setPreferredSize(new Dimension(100, 20));
        removeBtn.setMargin(new Insets(0, 0, 0, 0));
        removeBtn.setText("< Remove");
        removeBtn.addActionListener(this);
        addAllBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addAllBtn.setMaximumSize(new Dimension(100, 20));
        addAllBtn.setMinimumSize(new Dimension(100, 20));
        addAllBtn.setPreferredSize(new Dimension(100, 20));
        addAllBtn.setMargin(new Insets(0, 0, 0, 0));
        addAllBtn.setText("Add All >>");
        addAllBtn.addActionListener(this);
        transferPanel.setLayout(gridBagLayout1);
        selectedVolsList.setModel(selectedVolsListModel);
        availableVolsList.setModel(availableVolsListModel);
        selectedVolsList.addListSelectionListener(this);
        availableVolsList.addListSelectionListener(this);
        addBtn.setText("Add >");
        addBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addBtn.setMaximumSize(new Dimension(100, 20));
        addBtn.setMinimumSize(new Dimension(100, 20));
        addBtn.setPreferredSize(new Dimension(100, 20));
        addBtn.setMargin(new Insets(0, 0, 0, 0));
        addBtn.addActionListener(this);
        availableVolsScrollPane.setPreferredSize(new Dimension(150, 110));
        selectedVolsScrollPane.setPreferredSize(new Dimension(150, 110));

        ImageIcon icon = StsIcon.createIcon("upArrow.gif");
        upBtn.setIcon(icon);
        icon = StsIcon.createIcon("downArrow.gif");
        downBtn.setIcon(icon);

        upBtn.addActionListener(this);
        downBtn.addActionListener(this);

//        autoChk = new StsCheckbox();
//        autoChk.setText("Auto Scale");
//        autoChk.setSelected(false);
//        autoChk.setEnabled(true);
//        autoChk.addActionListener(this);

        histPanel.setLayout(gridBagLayout1);
        histLabel.setText("Null Limit Based on Selected PostStack3d");
        histLabel.setHorizontalAlignment(SwingConstants.CENTER);
        histLabel.setFont(new java.awt.Font("Dialog", 3, 11));
        minLabel.setFont(new java.awt.Font("Dialog", 0, 11));
        minLabel.setText("Minimum:");
        maxLabel.setFont(new java.awt.Font("Dialog", 0, 11));
        maxLabel.setText("Maximum:");
        histogramMinLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMinLabel.setForeground(Color.gray);
        histogramMinPctLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMinPctLabel.setForeground(Color.gray);
        histogramMaxLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMaxLabel.setForeground(Color.gray);
        histogramMaxPctLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMaxPctLabel.setForeground(Color.gray);

        histPanel.setBorder(BorderFactory.createEtchedBorder());
        this.add(selectionPanel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 4, 0, 5), 0, 0));
        transferPanel.add(removeBtn,        new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 0));
        transferPanel.add(addAllBtn,      new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 0));
        transferPanel.add(addBtn,      new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 0, 1), 0, 0));
        transferPanel.add(removeAllBtn,         new GridBagConstraints(0, 3, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 0));

        this.add(histPanel,  new GridBagConstraints(0, 1, 5, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 4, 4, 5), 0, 0));
        histPanel.add(histLabel, new GridBagConstraints(0, 0, 5, 1, 0.0, 0.0
           ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
//        histPanel.add(autoChk,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
//           ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 2, 2, 2), 0, 0));
        histPanel.add(minLabel,  new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
           ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 2, 2, 2), 0, 0));
        histPanel.add(dataMinBean,  new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
           ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(15, 2, 10, 0), 0, 0));
        histPanel.add(maxLabel,  new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0
           ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 3, 0, 0), 0, 0));
        histPanel.add(dataMaxBean,  new GridBagConstraints(4, 1, 1, 1, 1.0, 0.0
           ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(10, 3, 0, 0), 0, 0));
        histPanel.add(dataHistogram,  new GridBagConstraints(0, 2, 5, 2, 1.0, 1.0
           ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 0, 0, 0), 0, 0));

        selectionPanel.add(downBtn,  new GridBagConstraints(4, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 2, 5), 0, 0));
        selectionPanel.add(upBtn,  new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 2, 0), 0, 0));

        selectionPanel.add(availableVolsScrollPane,  new GridBagConstraints(0, 0, 2, 2, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 10, 0), 100, 0));
        selectionPanel.add(selectedVolsScrollPane,  new GridBagConstraints(3, 0, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 0, 0, 5), 100, 0));

        selectedVolsScrollPane.getViewport().add(selectedVolsList, null);
        availableVolsScrollPane.getViewport().add(availableVolsList, null);

        selectionPanel.add(transferPanel,  new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 0, 0, 0), 74, 74));

    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        int[] selectedIndices = null;

        Object source = e.getSource();

        selectedVolsList.removeListSelectionListener(this);

        // Add selected volumes
        if(source == addBtn)
        {
            selectedIndices = availableVolsList.getSelectedIndices();
            for(i=0; i< selectedIndices.length; i++)
            {
                if(selectedVolsListModel.indexOf(availableVolsList.getModel().getElementAt(selectedIndices[i])) >= 0)
                    continue;
                selectedVolsListModel.addElement(availableVolsList.getModel().getElementAt(selectedIndices[i]));
            }
        }
        // Remove selected volume(s)
        else if (source == removeBtn)
        {
            selectedIndices = selectedVolsList.getSelectedIndices();
            for(i= selectedIndices.length - 1; i>=0; i--)
                selectedVolsListModel.removeElementAt(selectedIndices[i]);
        }
        // Remove all volumes
        else if (source == removeAllBtn)
            selectedVolsListModel.removeAllElements();

        // Add all volumes
        else if (source == addAllBtn)
        {
            selectedVolsListModel.removeAllElements();
            for(i = 0; i < availableVolsList.getModel().getSize(); i++)
                selectedVolsListModel.addElement(availableVolsList.getModel().getElementAt(i));
        }
        // Up Button
        else if (source == upBtn)
        {
            if(selectedVol == null)
                new StsMessage(model.win3d, StsMessage.WARNING,"Must select item to move in left list.");
            else
            {
                moveSelectedItem(UP);
            }
            ;
        }
        // Down Button
        else if (source == downBtn)
        {
            if(selectedVol == null)
                new StsMessage(model.win3d, StsMessage.WARNING, "Must select item to move in left list.");
            else
            {
                moveSelectedItem(DOWN);
            }
        }
//        else if(source == autoChk)
//        {
//            setAutoChk(autoChk.isSelected());
//        }
        selectedVolsList.addListSelectionListener(this);
    }

    private void moveSelectedItem(int dir)
    {
        StsSeismicVolume[] list = null;
        int i;

        list = getVolumeList();
        if(list == null)
            return;
        switch(dir)
        {
            case UP:
                if(selectedIdx != 0)
                {
                    selectedVolsListModel.removeAllElements();
                    for (i = 0; i < selectedIdx-1; i++)
                        selectedVolsListModel.addElement(list[i].getName());
                    selectedVolsListModel.addElement(list[selectedIdx].getName());
                    selectedVolsListModel.addElement(list[selectedIdx-1].getName());
                    for (i = selectedIdx+1; i < list.length; i++)
                        selectedVolsListModel.addElement(list[i].getName());
                    setSelected(selectedIdx);
                }
                break;
            case DOWN:
                if(selectedIdx != (list.length-1))
                {
                    selectedVolsListModel.removeAllElements();
                    for (i = 0; i < selectedIdx; i++)
                        selectedVolsListModel.addElement(list[i].getName());
                    selectedVolsListModel.addElement(list[selectedIdx+1].getName());
                    selectedVolsListModel.addElement(list[selectedIdx].getName());
                    for (i = selectedIdx+2; i < list.length; i++)
                        selectedVolsListModel.addElement(list[i].getName());
                    setSelected(selectedIdx + 1);
                }

                break;
        }
    }

    private void setSelected(int idx)
    {
        selectedIdx = idx;
        selectedVol = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, (String) selectedVolsList.getModel().getElementAt(selectedIdx));
        selectedVolsList.setSelectedIndex(selectedIdx);
        resetLimits();
    }

    public void valueChanged(ListSelectionEvent e)
    {
        Object source = e.getSource();
        int[] selectedIndices = null;
        String volName = null;

        if(source == selectedVolsList)
        {
            selectedVolsList.removeListSelectionListener(this);
            selectedIndices = selectedVolsList.getSelectedIndices();
            if(selectedIndices != null)
                volName = (String) selectedVolsList.getModel().getElementAt(selectedIndices[selectedIndices.length - 1]);
            else
                volName = null;
            setSelected(selectedIndices[selectedIndices.length-1]);
            selectedVolsList.addListSelectionListener(this);
        }
        else if(source == availableVolsList)
            selectedVol = null;

        return;
    }

    public StsSeismicVolume[] getVolumeList()
    {
        String setName = null;
        StsSeismicVolume[] volSet = new StsSeismicVolume[selectedVolsList.getModel().getSize()];
        for(int i=0; i<selectedVolsList.getModel().getSize(); i++)
        {
            setName = (String)selectedVolsList.getModel().getElementAt(i);
            volSet[i] = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, setName);
        }
        return volSet;
    }

    private void resetLimits()
    {
        if(selectedVol == null)
            return;
//        if(autoChk.isSelected())
//            return;

        if((dataMin < selectedVol.getDataMin()) || (dataMin == StsParameters.largeFloat))
            dataMin = selectedVol.getDataMin();
        if((dataMax > selectedVol.getDataMax()) || (dataMax == StsParameters.smallFloat))
            dataMax = selectedVol.getDataMax();

        clearHistogramPanel();
        displayHistogram(selectedVol);

        dataMinBean.setValueAndRange(dataMin, selectedVol.getDataMin(), selectedVol.getDataMax());
        dataMaxBean.setValueAndRange(dataMax, selectedVol.getDataMin(), selectedVol.getDataMax());

        setDataMin(dataMin);
        setDataMax(dataMax);
    }

    public void displayHistogram(StsSeismicVolume volume)
    {
        histogram.updateData(volume.dataHist, volume.getDataMin(), volume.getDataMax());

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
    public float getDataMax() { return dataMax; }
    public void setDataMax(float value)
    {
        dataMax = value;
//        if(autoChk.isSelected())
//            return;
        if(histogram != null)
        {
            int minIdx = (int)((dataMin - selectedVol.getDataMin()) * (histogram.getNumberIndices()/(selectedVol.getDataMax() - selectedVol.getDataMin())));
            int maxIdx = (int)((dataMax - selectedVol.getDataMin()) * (histogram.getNumberIndices()/(selectedVol.getDataMax() - selectedVol.getDataMin())));
            histogram.setClip(minIdx, maxIdx);
            histogram.repaint();
            histogramMaxPctLabel.setText(new String(numberFormat.format(histogram.getTopPercentageClipped()) + "%"));
        }

    }
    public float getDataMin() { return dataMin; }
    public void setDataMin(float value)
    {
        dataMin = value;
//        if(autoChk.isSelected())
//            return;
        if(histogram != null)
        {
            int minIdx = (int)((dataMin - selectedVol.getDataMin()) * (histogram.getNumberIndices()/(selectedVol.getDataMax() - selectedVol.getDataMin())));
            int maxIdx = (int)((dataMax - selectedVol.getDataMin()) * (histogram.getNumberIndices()/(selectedVol.getDataMax() - selectedVol.getDataMin())));
            histogram.setClip(minIdx, maxIdx);
            histogram.repaint();
            histogramMinPctLabel.setText(new String(numberFormat.format(histogram.getBottomPercentageClipped()) + "%"));
        }

    }
/*
    public boolean getAutoChk() { return autoChk.isSelected(); }
    public void setAutoChk(boolean val)
    {
        if(val)
        {
            float[] range = ((StsCombinationVolumeWizard) wizard).setupVolume.calcDataRange();
            clearHistogramPanel();

            dataMinBean.setValueAndRange(range[0], range[0], range[1]);
            dataMaxBean.setValueAndRange(range[1], range[0], range[1]);

            setDataMin(range[0]);
            setDataMax(range[1]);
        }
        else
        {
            dataMin = StsParameters.largeFloat;
            dataMax = StsParameters.smallFloat;
            resetLimits();
        }
    }
*/
}
