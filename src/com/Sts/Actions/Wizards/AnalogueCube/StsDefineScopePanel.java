package com.Sts.Actions.Wizards.AnalogueCube;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineScopePanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsDefineScope wizardStep;
    private StsSubVolume[] availableSubVolumes = null;
    private StsSubVolume subVolume;
    private StsModel model = null;

    StsSeismicVolume sourceVolume = null;
    StsBoxSetSubVolume sourceBoxSet = null;
    StsBoxSubVolume sourceBox = null;

    private int nCurrentBoxRows = 10;
    private int nCurrentBoxCols = 10;
    private int nCurrentBoxSlices = 10;

    private float rowNumMin;
    private float rowNumMax;
    private float colNumMin;
    private float colNumMax;
    private float sliceMin;
    private float sliceMax;

    private boolean centerPicked = false;
    private StsCursorPoint cubePoint = null;

    ButtonGroup centerGroup = new ButtonGroup();
    StsJPanel beanPanel = StsJPanel.addInsets();
    private StsFloatFieldBean rowNumMinBean;
    private StsFloatFieldBean rowNumMaxBean;
    private JLabel rowLabel = new JLabel("Inline:");
    private StsFloatFieldBean colNumMinBean;
    private StsFloatFieldBean colNumMaxBean;
    private JLabel colLabel = new JLabel("Crossline:");
    private StsFloatFieldBean sliceMinBean;
    private StsFloatFieldBean sliceMaxBean;
    private JLabel zLabel = new JLabel("Z:");
    private StsGroupBox manualGroupBox = new StsGroupBox();

    private StsGroupBox centerGroupBox = new StsGroupBox();
    private JRadioButton maxBtn = new JRadioButton("Maximum");
    private JRadioButton minBtn = new JRadioButton("Minimum");
    private JRadioButton zeroBtn = new JRadioButton("Zero Crossing");
    private JRadioButton noneBtn = new JRadioButton("None");

    private GridBagLayout gridBagLayout = new GridBagLayout();
    private GridBagLayout gridBagLayout1 = new GridBagLayout();
    private GridBagLayout gridBagLayout2 = new GridBagLayout();
    JTextArea jTextArea1 = new JTextArea();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    public StsDefineScopePanel(StsWizard wizard, StsDefineScope defineScope)
    {
        this.wizard = wizard;
        this.wizardStep = defineScope;
        try
        {
            rowNumMinBean = new StsFloatFieldBean(this, "rowNumMin");
            rowNumMaxBean = new StsFloatFieldBean(this, "rowNumMax");
            colNumMinBean = new StsFloatFieldBean(this, "colNumMin");
            colNumMaxBean = new StsFloatFieldBean(this, "colNumMax");
            sliceMinBean = new StsFloatFieldBean(this, "sliceMin");
            sliceMaxBean = new StsFloatFieldBean(this, "sliceMax");

            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public boolean initialize(StsSeismicVolume sourceVolume)
    {
        model = wizard.getModel();
        this.sourceVolume = sourceVolume;
        model.setCurrentObjectDisplayAndToolbar(sourceVolume);

        if(sourceBoxSet == null)
        {
            sourceBoxSet = new StsBoxSetSubVolume(false);
			sourceBoxSet.setIsVisible(true);
            defineBox(null);
            sourceBoxSet.setActionEdit();
			model.addDisplayableInstance(sourceBoxSet);
        }
        else
        {
//            model.setCurrentObjectDisplayAndToolbar(getSourceBox().getVolume());
//            model.setCurrentObject(getSourceBox().getVolume());
//            model.win3d.getCursor3d().clearTextureDisplays();
//            model.win3dDisplay();
        }
        model.win3dDisplay();
        defineSourceVolume();
        return true;
    }

    public boolean defineSourceVolume()
    {
//        sourceVolume = (StsSeismicVolume) model.getCurrentObject(StsSeismicVolume.class);
        if (sourceVolume != null)
        {
            if(sourceBoxSet != null)
            {
                if(sourceBox != null)
                    sourceBox.setVolume(sourceVolume);
            }
            return true;
        }
        new StsMessage(model.win3d, StsMessage.WARNING, "No seismic volume currently defined as source volume.");
        return false;
    }

    void jbInit() throws Exception
    {
        manualGroupBox.setLabel("Manual Input");
        manualGroupBox.setForeground(Color.gray);
        manualGroupBox.setFont(new java.awt.Font("Dialog", 0, 11));
        manualGroupBox.setLayout(gridBagLayout);

        rowNumMinBean.setValue(0.0f);
        rowNumMinBean.setPreferredSize(new Dimension(100, 20));
        rowNumMinBean.setMinimumSize(new Dimension(40, 20));
        rowNumMinBean.setMaximumSize(new Dimension(150, 20));

        rowNumMaxBean.setValue(100.0f);
        rowNumMaxBean.setPreferredSize(new Dimension(100, 20));
        rowNumMaxBean.setMinimumSize(new Dimension(40, 20));
        rowNumMaxBean.setMaximumSize(new Dimension(150, 20));

        colNumMinBean.setValue(0.0f);
        colNumMinBean.setPreferredSize(new Dimension(100, 20));
        colNumMinBean.setMinimumSize(new Dimension(40, 20));
        colNumMinBean.setMaximumSize(new Dimension(150, 20));

        colNumMaxBean.setValue(100.0f);
        colNumMaxBean.setPreferredSize(new Dimension(100, 20));
        colNumMaxBean.setMinimumSize(new Dimension(40, 20));
        colNumMaxBean.setMaximumSize(new Dimension(150, 20));

        sliceMinBean.setValue(0.0f);
        sliceMinBean.setMinimum(0.0f);
        sliceMinBean.setPreferredSize(new Dimension(100, 20));
        sliceMinBean.setMinimumSize(new Dimension(40, 20));
        sliceMinBean.setMaximumSize(new Dimension(150, 20));

        sliceMaxBean.setValue(100.0f);
        sliceMaxBean.setPreferredSize(new Dimension(100, 20));
        sliceMaxBean.setMinimumSize(new Dimension(40, 20));
        sliceMaxBean.setMaximumSize(new Dimension(150, 20));

        centerGroupBox.setLabel("Center Box On ...");
        centerGroupBox.setForeground(Color.gray);
        centerGroupBox.setFont(new java.awt.Font("Dialog", 0, 11));
        centerGroupBox.setLayout(gridBagLayout2);

        maxBtn.addActionListener(this);
        minBtn.addActionListener(this);
        noneBtn.addActionListener(this);
        zeroBtn.addActionListener(this);

        zeroBtn.setSelected(true);

        this.setLayout(gridBagLayout1);
        beanPanel.setLayout(gridBagLayout3);
        jTextArea1.setBackground(Color.lightGray);
        jTextArea1.setFont(new java.awt.Font("Dialog", 1, 12));
        jTextArea1.setForeground(Color.blue);
        jTextArea1.setEditable(false);
        jTextArea1.setText("Define source subVolume graphically with the mouse and/or manually above. ");
        jTextArea1.setLineWrap(true);
        jTextArea1.setWrapStyleWord(true);
        this.add(beanPanel,  new GridBagConstraints(0, 0, 2, 3, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
        beanPanel.setBorder(BorderFactory.createEtchedBorder());
//        beanPanel.initializeLayout();

//        beanPanel.add(centerGroupBox, new GridBagConstraints(0, 1, 2, 1, 1.0, 1.0
//            , GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL, new Insets(2, 5, 0, 0), 0, 0));
        centerGroupBox.add(maxBtn, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 5, 0, 0), 0, 0));
        centerGroupBox.add(minBtn, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 5, 0, 0), 0, 0));
        centerGroupBox.add(zeroBtn, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 5, 0, 0), 0, 0));
        centerGroupBox.add(noneBtn, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 5, 0, 0), 0, 0));

        beanPanel.add(manualGroupBox,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 4), 288, 51));
        manualGroupBox.add(rowLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 5, 0, 0), 0, 0));
        manualGroupBox.add(rowNumMinBean, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 5, 0, 0), 0, 0));
        manualGroupBox.add(rowNumMaxBean, new GridBagConstraints(2, 0, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 5, 0, 0), 0, 0));
        manualGroupBox.add(colLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 4, 0, 0), 0, 0));
        manualGroupBox.add(colNumMinBean, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 4, 0, 0), 0, 0));
        manualGroupBox.add(colNumMaxBean, new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 4, 0, 0), 0, 0));
        manualGroupBox.add(zLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 4, 0, 0), 0, 0));
        manualGroupBox.add(sliceMinBean, new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 4, 0, 0), 0, 0));
        manualGroupBox.add(sliceMaxBean, new GridBagConstraints(2, 2, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 4, 0, 0), 0, 0));
        beanPanel.add(jTextArea1,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 48, 0, 56), 0, 11));

        centerGroup.add(zeroBtn);
        centerGroup.add(maxBtn);
        centerGroup.add(minBtn);
        centerGroup.add(noneBtn);
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        if (!centerPicked)
        {
            return defineCube(mouse, glPanel3d);
        }
        else
        {
            return setEditMode(mouse, glPanel3d);
        }
    }

    private boolean setEditMode(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        if(!editCube(mouse, glPanel3d))
        {
            StsCursor3d cursor3d = model.getGlPanel3d().getCursor3d();
            cubePoint = cursor3d.getCursorPoint(model.getGlPanel3d(), mouse);

            boolean createNewBox = StsYesNoDialog.questionValue(model.win3d, "Do you wish to make a new source box at this picked location?");
            if(createNewBox)
            {
                centerPicked = false;
                defineBox(cubePoint);
                return true;
            }
            else
            {
                sourceBoxSet.setCurrentBox(sourceBox);
                return true;
            }
        }
        else
            return true;

    }
    private boolean defineCube(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        if (mouse.getCurrentButton() != StsMouse.LEFT)
        {
            return false;
        }

        int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        // set temporary cursorPoint on mouse PRESS or DRAG
        if (buttonState == StsMouse.PRESSED || buttonState == StsMouse.DRAGGED)
        {
            StsView currentView = glPanel3d.getView();
            if (currentView instanceof StsView3d)
            {
                StsCursor3d cursor3d = model.getGlPanel3d().getCursor3d();
                cubePoint = cursor3d.getCursorPoint(model.getGlPanel3d(), mouse);
                if (cubePoint == null)
                {
                    return false;
                }
                cursor3d.setCurrentDirNo(cubePoint.dirNo);
            }
            else if (currentView instanceof StsViewCursor)
            {
                cubePoint = ( (StsViewCursor) currentView).getCursorPoint(mouse);
                if (cubePoint == null)
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
        // permanently add this point when mouse.RELEASED
        else if (buttonState == StsMouse.RELEASED)
        {
            defineSourceVolume();
            if (sourceVolume == null)
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "No seismic volumes loaded.");
                return false;
            }
            if (cubePoint == null)
            {
                return false;
            }
            defineBox(cubePoint);
        }
        model.win3dDisplay();
        return true;
    }

    private boolean editCube(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        boolean edited = sourceBoxSet.edit(mouse, glPanel3d);
        if (edited)
        {
            model.subVolumeChanged();
            updateBeansFromBox();
        }

        return edited;
    }

    private void updateBeansFromBox()
    {
        rowNumMinBean.setValue(sourceBox.rowNumMin);
        rowNumMaxBean.setValue(sourceBox.rowNumMax);
        colNumMinBean.setValue(sourceBox.colNumMin);
        colNumMaxBean.setValue(sourceBox.colNumMax);
        sliceMinBean.setValue(sourceBox.zMin);
        sliceMaxBean.setValue(sourceBox.zMax);
//        sliceMinBean.setValue(sourceBox.getSliceMin());
//        sliceMaxBean.setValue(sourceBox.getSliceMax());
    }

    public void actionPerformed(ActionEvent e)
    {
        if (noneBtn.isSelected())
        {
            ;
        }
        else if (maxBtn.isSelected())
        {
            ;
        }
        else if (minBtn.isSelected())
        {
            ;
        }
        else if (zeroBtn.isSelected())
        {
            ;
        }

    }

    public StsBoxSubVolume getSourceBox()
    {
		return sourceBox;
    }

	public StsBoxSetSubVolume getSourceBoxSet()
	{
		return sourceBoxSet;
	}

    public float getRowNumMin()
    {
        if (!boxOK())
        {
            return 0.0f;
        }
        return sourceBox.getRowNumMin();
    }

    public float getRowNumMax()
    {
        if (!boxOK())
        {
            return 0.0f;
        }
        return sourceBox.getRowNumMax();
    }

    public float getColNumMin()
    {
        if (!boxOK())
        {
            return 0.0f;
        }
        return sourceBox.getColNumMin();
    }

    public float getColNumMax()
    {
        if (!boxOK())
        {
            return 0.0f;
        }
        return sourceBox.getColNumMax();
    }

    public float getSliceMin()
    {
        if (!boxOK())
        {
            return 0.0f;
        }
        return sourceBox.getSliceMin();
    }

    public float getSliceMax()
    {
        if (!boxOK())
        {
            return 0.0f;
        }
        return sourceBox.getSliceMax();
    }

    private boolean boxOK()
    {
        return sourceBox != null;
    }

    public void setRowNumMin(float value)
    {
        sourceBox.setRowNumMin(value);
        rowNumMinBean.setValue(value);
        model.win3dDisplay();
    }

    public void setRowNumMax(float value)
    {
        sourceBox.setRowNumMax(value);
        rowNumMaxBean.setValue(value);
        model.win3dDisplay();
    }

    public void setColNumMin(float value)
    {
        sourceBox.setColNumMin(value);
        colNumMinBean.setValue(value);
        model.win3dDisplay();
   }

    public void setColNumMax(float value)
    {
        sourceBox.setColNumMax(value);
        colNumMaxBean.setValue(value);
        model.win3dDisplay();
    }

    public void setSliceMin(float value)
    {
        sliceMin = value;
        this.sliceMinBean.setValue(value);
        sourceBox.setSliceMin(sourceVolume.getSliceCoor(value));
        model.win3dDisplay();
   }

    public void setSliceMax(float value)
    {
        sliceMax = value;
        this.sliceMaxBean.setValue(value);
        sourceBox.setSliceMax(sourceVolume.getSliceCoor(value));
        model.win3dDisplay();
   }

    public void defineBox(StsCursorPoint centerPoint)
    {
        defineSourceVolume();
        nCurrentBoxRows = 10;
        nCurrentBoxCols = 10;
        nCurrentBoxSlices = 10;

        StsBoxSetSubVolumeClass boxSetClass = (StsBoxSetSubVolumeClass) model.getStsClass(StsBoxSetSubVolume.class);
        sourceBoxSet.deleteBoxes();
        boxSetClass.setCurrentObject(sourceBoxSet);

        if (centerPoint == null)
        {
            centerPicked = false;
            float rowNumCenter = (sourceVolume.rowNumMin + sourceVolume.rowNumMax) / 2;
            float colNumCenter = (sourceVolume.colNumMin + sourceVolume.colNumMax) / 2;
            rowNumMin = rowNumCenter - nCurrentBoxRows * sourceVolume.rowNumInc;
            rowNumMax = rowNumCenter + nCurrentBoxRows * sourceVolume.rowNumInc;
            colNumMin = colNumCenter - nCurrentBoxCols * sourceVolume.colNumInc;
            colNumMax = colNumCenter + nCurrentBoxCols * sourceVolume.colNumInc;
            float sliceCenter = (sourceVolume.getSliceMax() + sourceVolume.getSliceMax()) / 2;
            sliceMin = sliceCenter - nCurrentBoxSlices / 2;
            sliceMax = sliceCenter + nCurrentBoxSlices / 2;
            sourceBox = new StsBoxSubVolume(sourceVolume, rowNumMin, rowNumMax, colNumMin, colNumMax, sliceMin,
                                            sliceMax);
//            boxSetClass.setIsVisible(false);
        }
        else
        {
            centerPicked = true;
            sourceBox = new StsBoxSubVolume(centerPoint, sourceVolume, nCurrentBoxRows, nCurrentBoxCols,
                                            nCurrentBoxSlices);
//            boxSetClass.setIsVisible(true);
            model.subVolumeChanged(); // ??

        }

        rowNumMinBean.setValueAndRangeFixStep(sourceBox.rowNumMin, sourceVolume.rowNumMin, sourceVolume.rowNumMax,
                                              sourceVolume.getRowNumInc());
        rowNumMaxBean.setValueAndRangeFixStep(sourceBox.rowNumMax, sourceVolume.rowNumMin, sourceVolume.rowNumMax,
                                              sourceVolume.getRowNumInc());
        colNumMinBean.setValueAndRangeFixStep(sourceBox.colNumMin, sourceVolume.colNumMin, sourceVolume.colNumMax,
                                              sourceVolume.getColNumInc());
        colNumMaxBean.setValueAndRangeFixStep(sourceBox.colNumMax, sourceVolume.colNumMin, sourceVolume.colNumMax,
                                              sourceVolume.getColNumInc());
        sliceMinBean.setValueAndRangeFixStep(sourceBox.zMin, sourceVolume.zMin, sourceVolume.zMax, sourceVolume.zInc);
        sliceMaxBean.setValueAndRangeFixStep(sourceBox.zMax, sourceVolume.zMin, sourceVolume.zMax, sourceVolume.zInc);
//        sliceMinBean.setValueAndRangeFixStep(sourceBox.getSliceMin(), sourceVolume.getSliceMin(), sourceVolume.getSliceMax(), 1);
//        sliceMaxBean.setValueAndRangeFixStep(sourceBox.getSliceMax(), sourceVolume.getSliceMin(), sourceVolume.getSliceMax(), 1);

        sourceBoxSet.setStsColor(StsColor.CYAN);
        sourceBoxSet.add(sourceBox);
    }

    public void deleteBoxes()
    {
        if(sourceBoxSet == null) return;
        sourceBoxSet.deleteBoxes();
		model.removeDisplayableInstance(sourceBoxSet);
    }
}
