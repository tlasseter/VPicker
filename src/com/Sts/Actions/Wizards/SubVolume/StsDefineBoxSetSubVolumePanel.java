package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>e
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsDefineBoxSetSubVolumePanel extends StsFieldBeanPanel
{
    private StsWizard wizard;
    private StsDefineBoxSetSubVolume defineBoxSetSubVolume;

	private int nBoxRows;
	private int nBoxCols;
	private int nBoxSlices;

    StsJPanel beanPanel = new StsJPanel();
    StsGroupBox setPropertiesBox = new StsGroupBox("Box Set Properties");
    StsGroupBox propertiesBox = new StsGroupBox("Box Properties");
    StsGroupBox operationsBox = new StsGroupBox("Operations");
    public StsToggleButton pickButton;
    public StsToggleButton deleteButton;
    public StsToggleButton moveButton;
    public StsToggleButton editButton;
    static final String pickButtonTip = "Set to hexahedral center selection.";
    static final String deleteButtonTip = "Set to delete hexahedral.";
    static final String moveButtonTip = "Move the selected hexahedral.";
    static final String editButtonTip = "Stretch/shrink side of the selected hexahedral.";

    JLabel corRangeLabel = new JLabel(noCorrelString);
    JLabel waveLengthRangeLabel = new JLabel(noWaveLengthString);
    static final String noCorrelString = "Average Amplitude: not available.";
    static final String noWaveLengthString = "WaveLength Range: not available.";
    ButtonGroup buttonGroup = new ButtonGroup();

    StsColorComboBoxFieldBean cubeColorBean = new StsColorComboBoxFieldBean();
    StsIntFieldBean nBoxRowsBean = new StsIntFieldBean();
    StsIntFieldBean nBoxColsBean = new StsIntFieldBean();
    StsIntFieldBean nBoxSlicesBean = new StsIntFieldBean();

    StsFloatFieldBean xMinBean = new StsFloatFieldBean();
    StsFloatFieldBean xMaxBean = new StsFloatFieldBean();
    StsFloatFieldBean yMinBean = new StsFloatFieldBean();
    StsFloatFieldBean yMaxBean = new StsFloatFieldBean();
    StsFloatFieldBean zMinBean = new StsFloatFieldBean();
    StsFloatFieldBean zMaxBean = new StsFloatFieldBean();

    public StsDefineBoxSetSubVolumePanel(StsWizard wizard, StsDefineBoxSetSubVolume defineBoxSetSubVolume)
    {
        this.wizard = wizard;
        this.defineBoxSetSubVolume = defineBoxSetSubVolume;
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    void jbInit() throws Exception
    {

        setLayout(new GridBagLayout());
        add(operationsBox, new GridBagConstraints(0, 0, 2, 1, 1.0, 0.0,
            GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 5, 5), 0, 0));
        add(setPropertiesBox, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0,
            GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 5, 5), 0, 0));
        add(propertiesBox, new GridBagConstraints(0, 2, 1, 1, 1.0, 0.0,
            GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 5, 5), 0, 0));

        cubeColorBean.initializeColors(defineBoxSetSubVolume, "stsColor", "Color:", StsColor.colors8);

        int rows = wizard.getModel().getProject().getRotatedBoundingBox().getNRows();
        int cols = wizard.getModel().getProject().getRotatedBoundingBox().getNCols();
        int slices = wizard.getModel().getProject().getRotatedBoundingBox().getNSlices();
	    nBoxRows = rows/20;
		nBoxCols = cols/20;
		nBoxSlices = slices/20;
        nBoxRowsBean.initialize(this, "nBoxRows", 1, rows, "Inlines");
        nBoxColsBean.initialize(this, "nBoxCols", 1, cols, "Crosslines");
        nBoxSlicesBean.initialize(this, "nBoxSlices", 1, slices, "Slices");

        setPropertiesBox.add(cubeColorBean);
        setPropertiesBox.add(nBoxRowsBean);
        setPropertiesBox.add(nBoxColsBean);
        setPropertiesBox.add(nBoxSlicesBean);

        xMinBean.initialize(defineBoxSetSubVolume, "boxXMin", true, "min X");
        xMaxBean.initialize(defineBoxSetSubVolume, "boxXMax", true, "max X");
        yMinBean.initialize(defineBoxSetSubVolume, "boxYMin", true, "min Y");
        yMaxBean.initialize(defineBoxSetSubVolume, "boxYMax", true, "max Y");
        zMinBean.initialize(defineBoxSetSubVolume, "boxZMin", true, "min Z");
        zMaxBean.initialize(defineBoxSetSubVolume, "boxZMax", true, "max Z");

        propertiesBox.add(xMinBean);
        propertiesBox.add(xMaxBean);
        propertiesBox.add(yMinBean);
        propertiesBox.add(yMaxBean);
        propertiesBox.add(zMinBean);
        propertiesBox.add(zMaxBean);

        pickButton = new StsToggleButton("Pick", pickButtonTip, defineBoxSetSubVolume, "setDefineCube");
        moveButton = new StsToggleButton("Move", moveButtonTip, defineBoxSetSubVolume, "setMoveCube");
        editButton = new StsToggleButton("Edit", editButtonTip, defineBoxSetSubVolume, "setEditCube");
        deleteButton = new StsToggleButton("Delete", deleteButtonTip, defineBoxSetSubVolume, "setDeleteCube");

        operationsBox.add(pickButton, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
        operationsBox.add(moveButton, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
        operationsBox.add(editButton, new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
        operationsBox.add(deleteButton, new GridBagConstraints(3, 0, 1, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));

        operationsBox.add(corRangeLabel, new GridBagConstraints(0, 1, 3, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 5, 5, 5), 0, 0));
        operationsBox.add(waveLengthRangeLabel, new GridBagConstraints(0, 2, 3, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 5, 5, 5), 0, 0));

        buttonGroup.add(pickButton);
        buttonGroup.add(moveButton);
        buttonGroup.add(editButton);
        buttonGroup.add(deleteButton);
    }

    public void initializeButtonGroup()
    {
        pickButton.setSelected(true);
    }

    public void setStsColor(StsColor color) { cubeColorBean.setStsColor(color); }
    public StsColor getStsColor() { return cubeColorBean.getStsColor(); }

    public void setNBoxRows(int nRows) { nBoxRows = nRows; }
    public void setNBoxCols(int nCols) { nBoxCols = nCols; }
    public void setNBoxSlices(int nSlices) { nBoxSlices = nSlices; }
    public int getNBoxRows() { return nBoxRows; }
    public int getNBoxCols() { return nBoxCols; }
    public int getNBoxSlices() { return nBoxSlices; }

    public void updateBoxProperties()
    {
        cubeColorBean.getValueFromPanelObject();
        xMinBean.getValueFromPanelObject();
        xMaxBean.getValueFromPanelObject();
        yMinBean.getValueFromPanelObject();
        yMaxBean.getValueFromPanelObject();
        zMinBean.getValueFromPanelObject();
        zMaxBean.getValueFromPanelObject();
    }
}
