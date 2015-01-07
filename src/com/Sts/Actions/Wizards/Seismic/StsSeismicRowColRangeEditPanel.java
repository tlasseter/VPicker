package com.Sts.Actions.Wizards.Seismic;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.Lasseter
 * @version 1.1
 */

public class StsSeismicRowColRangeEditPanel extends StsFieldBeanPanel
{
    EditedBoundingBox editedBoundingBox = new EditedBoundingBox();
    StsSeismicBoundingBox[] volumes;
    private boolean rowNumInverted = false;
    private StsGroupBox groupBox;
    private StsFloatFieldBean rowNumMinBean;
    private StsFloatFieldBean rowNumMaxBean;
    private StsFloatFieldBean rowNumIncBean;
    private StsFloatFieldBean colNumMinBean;
    private StsFloatFieldBean colNumMaxBean;
    private StsFloatFieldBean colNumIncBean;
	private JLabel minLabel = new JLabel("Minimum");
	private JLabel maxLabel = new JLabel("Maximum");
	private JLabel intervalLabel = new JLabel("Increment");

	public StsSeismicRowColRangeEditPanel(String boxTitle)
	{
        try
		{
            constructBeans();
            buildPanel(boxTitle);
        }
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

    /** fix edit ranges so cannot extrapolate or subdivide */
    public void initialize(StsSeismicBoundingBox[] volumes)
	{
        this.volumes = volumes;
        editedBoundingBox.initialize(volumes);
        editedBoundingBox.initializeBeanValuesAndRanges();
        setPanelObject(editedBoundingBox);
    }

    /** Allow row-col range to be extended and subdivided, i.e., no restrictions on bean range and step. */
    public void initializeExtend(StsSeismicBoundingBox[] volumes)
	{
        this.volumes = volumes;
        editedBoundingBox.initialize(volumes);
        setPanelObject(editedBoundingBox);
    }

    private void constructBeans()
    {
        rowNumMinBean = new StsFloatFieldBean(EditedBoundingBox.class, "rowNumMin", "Inline:");
        rowNumMaxBean = new StsFloatFieldBean(EditedBoundingBox.class, "rowNumMax");
        rowNumIncBean = new StsFloatFieldBean(EditedBoundingBox.class, "rowNumInc");
        colNumMinBean = new StsFloatFieldBean(EditedBoundingBox.class, "colNumMin", "Crossline:");
        colNumMaxBean = new StsFloatFieldBean(EditedBoundingBox.class, "colNumMax");
        colNumIncBean = new StsFloatFieldBean(EditedBoundingBox.class, "colNumInc");
    }

    private void buildPanel(String boxTitle)
	{
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.NORTH;

        groupBox = new StsGroupBox(boxTitle);
        groupBox.gbc.gridx = 1; // leave empty cell
        groupBox.addToRow(minLabel);
        groupBox.addToRow(maxLabel);
        groupBox.addEndRow(intervalLabel);

        groupBox.addToRow(rowNumMinBean);
        groupBox.addToRow(rowNumMaxBean);
        groupBox.addEndRow(rowNumIncBean);

        groupBox.addToRow(colNumMinBean);
        groupBox.addToRow(colNumMaxBean);
        groupBox.addEndRow(colNumIncBean);
        add(groupBox);

    }

    public StsEditedBoundingBox getEditedBox()
    {
        StsEditedBoundingBox[] editedBoundingBoxes = getEditedBoxes();
        if(editedBoundingBoxes == null) return null;
        return editedBoundingBoxes[0];
    }
    public StsEditedBoundingBox[] getEditedBoxes()
    {
        // create a crop box for each volume which describes the relationship between the original volume and the cropped range
        int nVolumes = volumes.length;
        StsEditedBoundingBox[] editedBoundingBoxes = new StsEditedBoundingBox[nVolumes];
        for(int i=0; i<volumes.length; i++)
            editedBoundingBoxes[i] = new StsEditedBoundingBox(volumes[i], editedBoundingBox, false);
        if(editedBoundingBoxes == null) return null;
        return editedBoundingBoxes;
    }

    class EditedBoundingBox extends StsEditedBoundingBox
    {
        EditedBoundingBox()
        {
            super(false);
        }

        private void initialize(StsSeismicBoundingBox[] volumes)
        {
            StsSeismicBoundingBox firstVolume = volumes[0];
            initializeOriginAndAngle(firstVolume);

            for(int n = 0; n < volumes.length; n++)
                addBoundingBox(volumes[n]);
        }

        private void initializeBeanValuesAndRanges()
        {
            rowNumInverted = rowNumMin > rowNumMax;
            if(!rowNumInverted)
            {
                rowNumMinBean.setValueAndRangeFixStep(rowNumMin, rowNumMin, rowNumMax, rowNumInc);
                rowNumMaxBean.setValueAndRangeFixStep(rowNumMax, rowNumMin, rowNumMax, rowNumInc);
                rowNumIncBean.setValueAndRangeFixStep(rowNumInc, rowNumInc, rowNumMax - rowNumMin, rowNumInc);

            }
            else
            {
                rowNumMinBean.setValueAndRangeFixStep(rowNumMax, rowNumMax, rowNumMin, -rowNumInc);
                rowNumMaxBean.setValueAndRangeFixStep(rowNumMin, rowNumMax, rowNumMin, -rowNumInc);
                rowNumIncBean.setValueAndRangeFixStep(-rowNumInc, -rowNumInc, rowNumMin - rowNumMax, -rowNumInc);
            }

            colNumMinBean.setValueAndRangeFixStep(colNumMin, colNumMin, colNumMax, colNumInc);
            colNumMaxBean.setValueAndRangeFixStep(colNumMax, colNumMin, colNumMax, colNumInc);
            colNumIncBean.setValueAndRangeFixStep(colNumInc, colNumInc, colNumMax - colNumMin, colNumInc);
        }

        public void setRowNumMin(float rowNumMin)
        {
            if(!rowNumInverted)
                this.rowNumMin = rowNumMin;
            else
                this.rowNumMax = rowNumMin;
        }
        public float getRowNumMin()
        {
            if(!rowNumInverted)
                return rowNumMin;
            else
                return rowNumMax;
        }

        public void setRowNumMax(float rowNumMax)
        {
            if(!rowNumInverted)
                this.rowNumMax = rowNumMax;
            else
                this.rowNumMin = rowNumMax;
        }

        public float getRowNumMax()
        {
            if(!rowNumInverted)
                return rowNumMax;
            else
                return rowNumMin;
        }

        public void setRowNumInc(float value)
        {
            rowNumMinBean.fixStep(value);
            rowNumMaxBean.fixStep(value);
            if(rowNumInverted) value = -value;
            rowNumInc = value;
            rowNumMax = Math.round((rowNumMax - rowNumMin) / rowNumInc) * rowNumInc + rowNumMin;
            if(!rowNumInverted)
                rowNumMaxBean.setValue(rowNumMax);
            else
                rowNumMinBean.setValue(rowNumMax);
            nRows = Math.round((rowNumMax - rowNumMin)/rowNumInc) + 1;
        }

        public float getRowNumInc()
        {
            if(!rowNumInverted)
                return rowNumInc;
            else
                return -rowNumInc;
        }

        public void setColNumInc(float value)
        {
            colNumInc = value;
            colNumMinBean.fixStep(colNumInc);
            colNumMaxBean.fixStep(colNumInc);
            colNumMax = Math.round((colNumMax - colNumMin) / colNumInc) * colNumInc + colNumMin;
            colNumMaxBean.setValue(colNumMax);
            nCols = Math.round((colNumMax - colNumMin)/colNumInc) + 1;
        }
    }

    static public void main(String[] args)
    {
        StsModel model = StsModel.constructor();
        StsSeismicRowColRangeEditPanel panel = new StsSeismicRowColRangeEditPanel("Box Title");
        StsToolkit.createDialog(panel);
	}
}