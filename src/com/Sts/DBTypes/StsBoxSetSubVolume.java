package com.Sts.DBTypes;

import com.Sts.Actions.Wizards.SubVolume.*;
import com.Sts.DB.DBCommand.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Reflect.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;
import java.text.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsBoxSetSubVolume extends StsSubVolume implements StsTreeObjectI, StsInstance3dDisplayable
{
	private StsObjectRefList boxes;
    protected StsColor stsColor = new StsColor(StsColor.RED);
    transient boolean isEditing = false;
	transient StsBoxSubVolume currentBox = null;
	transient int pickedPointIndex = -1;
	transient byte action = ACTION_NONE;
	transient StsSeismicVolume currentVolume = null;
	transient boolean honorVolume = false;
	transient boolean displayVoxels = false;

	static byte ACTION_NONE = StsDefineBoxSetSubVolume.ACTION_NONE;
	static byte ACTION_MOVE_BOX = StsDefineBoxSetSubVolume.ACTION_MOVE;
	static byte ACTION_MOVE_POINT = StsDefineBoxSetSubVolume.ACTION_EDIT;
	static byte ACTION_DELETE = StsDefineBoxSetSubVolume.ACTION_DELETE;

	static StsObjectPanel objectPanel = null;

    static public StsFieldBean[] displayFields = null;

	static public final StsFieldBean[] propertyFields = null;

	static StsJOGLPick pickCurrentBoxSetCenterPoint = null;
	static StsJOGLPick pickCurrentBoxSetFacePoint = null;

	public StsBoxSetSubVolume()
	{
	}

	public StsBoxSetSubVolume(boolean persistent)
	{
		super(persistent);
		initializeVisibleFlags();
	}

	public StsBoxSetSubVolume(String name)
	{
		this(true, name);
	}

	public StsBoxSetSubVolume(boolean persistent, String name)
	{
		super(persistent);
		setName(name);
		initializeVisibleFlags();
	}

	public void add(StsBoxSubVolume box)
	{
        if(boxes == null)
        {
            boxes = StsObjectRefList.constructor(4, 4, "boxes", this);
        }
        boxes.add(box);
		currentBox = box;
	}

	public void setDisplayVoxels(boolean displayVoxels)
	{
		this.displayVoxels = displayVoxels;
	}

	public boolean getDisplayVoxels() { return displayVoxels; }

	public boolean anyDependencies()
	{
		return false;
	}

	public StsFieldBean[] getDisplayFields()
	{
        if(displayFields != null) return displayFields;
        displayFields = new StsFieldBean[]
		{
            new StsStringFieldBean(StsBoxSetSubVolume.class, "name", "Name"),
            new StsBooleanFieldBean(StsBoxSetSubVolume.class, "isVisible", "Visible"),
    //		new StsBooleanFieldBean(StsBoxSetSubVolume.class, "displayVoxels", "Voxels"),
            new StsBooleanFieldBean(StsBoxSetSubVolume.class, "isApplied", "Applied"),
    //		new StsBooleanFieldBean(StsBoxSetSubVolume.class, "isInclusive", "Inclusive"),
            new StsColorComboBoxFieldBean(StsBoxSetSubVolume.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors())
        };
        return displayFields;
	}

    public StsFieldBean[] getPropertyFields()
    {
        return null;
    }

    public Object[] getChildren()
	{
        return new Object[0];
	}

	public StsObjectPanel getObjectPanel()
	{
		if (objectPanel == null)
		{
			objectPanel = StsObjectPanel.constructor(this, true);
		}
		return objectPanel;
	}

    public boolean delete()
    {
        boxes.deleteAll();
        boolean success = super.delete();
        currentModel.viewObjectChangedAndRepaint(this, this);
        return success;
    }

    public void display(StsGLPanel glPanel, boolean isCurrentObject)
	{
		if (!isVisible) return;

        if(boxes == null) return;
        int nBoxes = boxes.getSize();
		for (int n = 0; n < nBoxes; n++)
		{
			StsBoxSubVolume box = (StsBoxSubVolume) boxes.getElement(n);
			boolean editing = currentBox == box;
			if(displayVoxels)
				box.displayVoxels(glPanel);
			else
				box.display(glPanel, stsColor, action, editing);
		}
	}

	public void display(StsGLPanel3d glPanel3d)
	{
		display(glPanel3d, false);
	}

	public void setStsColor(StsColor stsColor)
	{
        this.stsColor = stsColor;
        currentModel.addTransactionCmd("SubVolume color change", new StsChangeCmd(this, stsColor, "stsColor", false));
        currentModel.win3dDisplayAll();
    }

	public StsColor getStsColor()
	{
        return stsColor;
	}

	public StsGridBoundingBox getGridBoundingBox()
	{
		StsGridBoundingBox gridBoundingBox = new StsGridBoundingBox(false);
		int nBoxes = boxes.getSize();
		for (int n = 0; n < nBoxes; n++)
		{
			StsBoxSubVolume box = (StsBoxSubVolume) boxes.getElement(n);
			gridBoundingBox.addBoundingBox(box.getGridBoundingBox());
		}
		return gridBoundingBox;
	}

	public void setCurrentVolume(StsSeismicVolume vol)
	{
		currentVolume = vol;
	}

	public void setHonorVolumes(boolean honor)
	{
		honorVolume = honor;
	}

	public void addUnion(byte[] subVolumePlane, int dir, float dirCoordinate, StsSeismicVolume vol, boolean honorVols, byte zDomainData)
	{
		currentVolume = vol;
		honorVolume = honorVols;
		addUnion(subVolumePlane, dir, dirCoordinate, vol, zDomainData);
	}

	public void addUnion(byte[] subVolumePlane, int dir, float dirCoordinate, StsRotatedGridBoundingBox cursor3dBoundingBox, byte zDomainData)
	{
		int nTotalCursorRows, nTotalCursorCols;

		if (!isApplied) return;

        if(boxes == null) return;

        int nBoxes = boxes.getSize();
		switch (dir)
		{
			case StsCursor3d.XDIR:
				nTotalCursorCols = cursor3dBoundingBox.getNSlices();
				for (int n = 0; n < nBoxes; n++)
				{
					StsBoxSubVolume box = (StsBoxSubVolume) boxes.getElement(n);
					if ( (box.getVolume() != currentVolume) && (honorVolume))
					{
						continue;
					}
					if (dirCoordinate < box.xMin || dirCoordinate > box.xMax)
					{
						continue;
					}
					int cursorRowStart = StsMath.ceiling(cursor3dBoundingBox.getBoundedRowCoor(box.yMin));
					int cursorRowEnd = StsMath.floor(cursor3dBoundingBox.getBoundedRowCoor(box.yMax));
					int cursorColStart = StsMath.ceiling(cursor3dBoundingBox.getBoundedSliceCoor(box.zMin));
					int cursorColEnd = StsMath.floor(cursor3dBoundingBox.getBoundedSliceCoor(box.zMax));

					for (int cursorRow = cursorRowStart; cursorRow <= cursorRowEnd; cursorRow++)
					{
						int cursorIndex = cursorRow * nTotalCursorCols + cursorColStart;
						for (int cursorCol = cursorColStart; cursorCol <= cursorColEnd; cursorCol++)
						{
							subVolumePlane[cursorIndex++] = (byte) 1;
						}
					}
				}
				break;
			case StsCursor3d.YDIR:
				nTotalCursorCols = cursor3dBoundingBox.getNSlices();
				for (int n = 0; n < nBoxes; n++)
				{
					StsBoxSubVolume box = (StsBoxSubVolume) boxes.getElement(n);
					if ( (box.getVolume() != currentVolume) && (honorVolume))
					{
						continue;
					}
					if (dirCoordinate < box.yMin || dirCoordinate > box.yMax)
					{
						continue;
					}
					int cursorRowStart = StsMath.ceiling(cursor3dBoundingBox.getBoundedColCoor(box.xMin));
					int cursorRowEnd = StsMath.floor(cursor3dBoundingBox.getBoundedColCoor(box.xMax));
					int cursorColStart = StsMath.ceiling(cursor3dBoundingBox.getBoundedSliceCoor(box.zMin));
					int cursorColEnd = StsMath.floor(cursor3dBoundingBox.getBoundedSliceCoor(box.zMax));

					for (int cursorRow = cursorRowStart; cursorRow <= cursorRowEnd; cursorRow++)
					{
						int cursorIndex = cursorRow * nTotalCursorCols + cursorColStart;
						for (int cursorCol = cursorColStart; cursorCol <= cursorColEnd; cursorCol++)
						{
							subVolumePlane[cursorIndex++] = (byte) 1;
						}
					}
				}
				break;
			case StsCursor3d.ZDIR:
				nTotalCursorCols = cursor3dBoundingBox.getNCols();
				for (int n = 0; n < nBoxes; n++)
				{
					StsBoxSubVolume box = (StsBoxSubVolume) boxes.getElement(n);
					if ( (box.getVolume() != currentVolume) && (honorVolume))
					{
						continue;
					}
					if (dirCoordinate < box.zMin || dirCoordinate > box.zMax)
					{
						continue;
					}
					int cursorRowStart = StsMath.ceiling(cursor3dBoundingBox.getBoundedRowCoor(box.yMin));
					int cursorRowEnd = StsMath.floor(cursor3dBoundingBox.getBoundedRowCoor(box.yMax));
					int cursorColStart = StsMath.ceiling(cursor3dBoundingBox.getBoundedColCoor(box.xMin));
					int cursorColEnd = StsMath.floor(cursor3dBoundingBox.getBoundedColCoor(box.xMax));

					for (int cursorRow = cursorRowStart; cursorRow <= cursorRowEnd; cursorRow++)
					{
						int cursorIndex = cursorRow * nTotalCursorCols + cursorColStart;
						for (int cursorCol = cursorColStart; cursorCol <= cursorColEnd; cursorCol++)
						{
							subVolumePlane[cursorIndex++] = (byte) 1;
						}
					}
				}
				break;
		}
	}

	public void setIsEditing(boolean isEditing)
	{this.isEditing = isEditing;
	}

	public int getNBoxes()
	{return boxes.getSize();
	}

	public void setAction(byte action)
	{
		this.action = action;
	}

	public void setActionEdit()
	{
		action = ACTION_MOVE_POINT;
	}

	public boolean move(StsMouse mouse, StsGLPanel3d glPanel3d)
	{
		if (mouse.getCurrentButton() != StsMouse.LEFT)
		{
			return false;
		}
		int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

		if (buttonState == StsMouse.CLEARED)
		{
			return true;
		}

		checkInitializePickMethods(glPanel3d);

		if (buttonState == StsMouse.PRESSED)
		{
			currentBox = pickCenter(glPanel3d);
			if (currentBox != null)
			{
				action = ACTION_MOVE_BOX;
			}
			else
			{
				action = ACTION_NONE;
			}
			currentModel.win3dDisplay();
		}
		else if (buttonState == StsMouse.DRAGGED)
		{
			if (currentBox == null)
			{
				return true;
			}
			if (action != ACTION_MOVE_BOX)
			{
				StsMessageFiles.errorMessage("StsBoxSetSubVolume.action not set to ACTION_MOVE_BOX.");
				return true;
			}
			currentBox.moveBox(mouse, glPanel3d);
			currentModel.win3dDisplay();
		}
		else if (buttonState == StsMouse.RELEASED)
		{
			action = ACTION_NONE;
		}
		return true;
	}

	private void checkInitializePickMethods(StsGLPanel3d glPanel3d)
	{
		if (pickCurrentBoxSetCenterPoint == null)
		{
			GL gl = glPanel3d.getGL();
			StsMethod method = new StsMethod(StsBoxSetSubVolume.class, "pickCurrentBoxSetCenterPoint", gl, GL.class);
			pickCurrentBoxSetCenterPoint = new StsJOGLPick(glPanel3d, method, 10, StsMethodPick.PICK_CLOSEST);
		}
		if (pickCurrentBoxSetFacePoint == null)
		{
			GL gl = glPanel3d.getGL();
			StsMethod method = new StsMethod(StsBoxSetSubVolume.class, "pickCurrentBoxSetFacePoint", gl, GL.class);
			pickCurrentBoxSetFacePoint = new StsJOGLPick(glPanel3d, method, 10, StsMethodPick.PICK_CLOSEST);
		}
	}

	public boolean edit(StsMouse mouse, StsGLPanel3d glPanel3d)
	{
		if (mouse.getCurrentButton() != StsMouse.LEFT)
		{
			return true;
		}
		int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

		if (buttonState == StsMouse.CLEARED)
		{
			return true;
		}

		checkInitializePickMethods(glPanel3d);

		if (buttonState == StsMouse.PRESSED)
		{
			if (action != ACTION_MOVE_POINT)
			{
				StsBoxSubVolume pickedBox = pickCenter(glPanel3d);
				if (pickedBox == null)
				{
					currentBox = null;
					action = ACTION_NONE;
				}
				else
				{
					currentBox = pickedBox;
					action = ACTION_MOVE_POINT;
				}
			}
			else // action == ACTION_MOVE_POINT
			{
				boolean facePicked = pickFace(glPanel3d);
				if (!facePicked)
				{
					StsBoxSubVolume pickedBox = pickCenter(glPanel3d);
					if (pickedBox != null && currentBox != pickedBox)
					{
						currentBox = pickedBox;
					}
					else if (pickedBox == null)
					{
						currentBox = null;
						action = ACTION_NONE;
					}
				}
			}
			currentModel.win3dDisplay();
			return currentBox != null;
		}
		else if (buttonState == StsMouse.DRAGGED)
		{
			if (currentBox == null)
			{
				return true;
			}
			if (action != ACTION_MOVE_POINT)
			{
				StsMessageFiles.errorMessage("StsBoxSetSubVolume.action not set to ACTION_MOVE_POINT.");
				return true;
			}
			currentBox.movePoint(mouse, glPanel3d, pickedPointIndex);
			currentModel.win3dDisplay();
		}
		else if (buttonState == StsMouse.RELEASED)
		{
            action = ACTION_NONE;
		}
		return true;
	}

	public boolean delete(StsMouse mouse, StsGLPanel3d glPanel3d)
	{
		if (mouse.getCurrentButton() != StsMouse.LEFT)
		{
			return false;
		}
		int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

		if (buttonState == StsMouse.CLEARED)
		{
			return false;
		}

		checkInitializePickMethods(glPanel3d);

		if (buttonState == StsMouse.RELEASED)
		{
			currentBox = pickCenter(glPanel3d);
			if (currentBox == null)
			{
				new StsMessage(currentModel.win3d, StsMessage.WARNING, "No hexahedral selected.  Try again.");
				return false;
			}
			else
			{
				boxes.delete(currentBox);
				currentBox.delete();
				currentModel.win3dDisplay();
				return true;
			}
		}
		else
		{
			return false;
		}
	}

	private StsBoxSubVolume pickCenter(StsGLPanel3d glPanel3d)
	{
		try
		{
			if (!pickCurrentBoxSetCenterPoint.methodPick3d())
			{
                pickedPointIndex = StsBoundingBox.NONE;
                return null;
			}
			StsPickItem pickItem = pickCurrentBoxSetCenterPoint.pickItems[0];
			int nBox = pickItem.names[0];
			return (StsBoxSubVolume) boxes.getElement(nBox);
		}
		catch (Exception e)
		{return null;
		}
	}

	private boolean pickFace(StsGLPanel3d glPanel3d)
	{
		try
		{
			if (!pickCurrentBoxSetFacePoint.methodPick3d())
			{
                pickedPointIndex = StsBoundingBox.NONE;
                return false;
			}
			StsPickItem pickItem = pickCurrentBoxSetFacePoint.pickItems[0];
			int nBox = pickItem.names[0];
			StsBoxSubVolume pickedBox = (StsBoxSubVolume) boxes.getElement(nBox);
			if (pickedBox != currentBox)
			{
				currentBox = pickedBox;
				return true;
			}
			pickedPointIndex = pickItem.names[1];
			currentModel.win3dDisplay();
			return true;
		}
		catch (Exception e)
		{
            return false;
		}
	}

	static public void pickCurrentBoxSetCenterPoint(GL gl)
	{
		StsBoxSetSubVolumeClass boxSetClass = (StsBoxSetSubVolumeClass) StsObject.getCurrentModel().getStsClass(StsBoxSetSubVolume.class);
		StsBoxSetSubVolume currentBoxSet = (StsBoxSetSubVolume) boxSetClass.getCurrentObject();
		if (currentBoxSet == null)
		{
			return;
		}
		currentBoxSet.pickCenterPoint(gl);
	}

	public void pickCenterPoint(GL gl)
	{
		int nBoxes = boxes.getSize();
		for (int n = 0; n < nBoxes; n++)
		{
			StsBoxSubVolume box = (StsBoxSubVolume) boxes.getElement(n);
			boolean editing = currentBox == box;
			box.pickCenterPoint(gl, editing, n);
		}
	}

	static public void pickCurrentBoxSetFacePoint(GL gl)
	{
		StsBoxSetSubVolumeClass boxSetClass = (StsBoxSetSubVolumeClass) StsObject.getCurrentModel().getStsClass(StsBoxSetSubVolume.class);
		StsBoxSetSubVolume currentBoxSet = (StsBoxSetSubVolume) boxSetClass.getCurrentObject();
		if (currentBoxSet == null)
		{
			return;
		}
		currentBoxSet.pickFacePoint(gl);
	}

	public void pickFacePoint(GL gl)
	{
		int nBoxes = boxes.getSize();
		for (int n = 0; n < nBoxes; n++)
		{
			StsBoxSubVolume box = (StsBoxSubVolume) boxes.getElement(n);
			boolean editing = currentBox == box;
			box.pickFacePoint(gl, n);
		}
	}

	public void deleteBoxes()
	{
		for (int i = 0; i < boxes.getSize(); i++)
		{
			StsBoxSubVolume box = (StsBoxSubVolume) boxes.getElement(i);
			boxes.delete(box);
			box.delete();
		}
	}

	public void setCurrentBox(StsBoxSubVolume box)
	{currentBox = box;
	}

	public StsBoxSubVolume getCurrentBox()
	{return currentBox;
	}

	public StsObjectRefList getBoxes()
	{return boxes;
	}

    public boolean canExport() { return true; }
	public boolean export()
	{
		StsProject project = currentModel.getProject();
		StsBoxSubVolume currentBox = null;
		StsGridPoint center = null;
		DecimalFormat floatFormat = new DecimalFormat("#.00");
		String filename = null;

		StsBoxSetExportDialog exportDialog = new StsBoxSetExportDialog(currentModel.win3d, "Box Set Export Utility", true, this);
		exportDialog.setVisible(true);
		if (exportDialog.getMode() == exportDialog.CANCELED)
		{
			return false;
		}
		if (exportDialog.getFormat() == exportDialog.ASCII)
		{
			filename = project.getRootDirString() + getName() + ".txt";
		}
		else
		{
			filename = project.getRootDirString() + getName() + ".xml";

		}
		File file = new File(filename);
		if (file.exists())
		{
			boolean overWrite = StsYesNoDialog.questionValue(currentModel.win3d,
				"File " + filename + " already exists. Do you wish to overwrite it?");
			if (!overWrite)
			{
				return false;
			}
		}

		// ASCII Output
		if (exportDialog.getFormat() == exportDialog.ASCII)
		{
			StsFile out = StsFile.constructor(filename);
			if (exportDialog.getScope() == exportDialog.CENTERS)
			{
				String[] lines = new String[getNBoxes() + 5];
				lines[0] = "---------------------------------------------------------";
				lines[1] = " Exported from S2S Systems ";
				lines[2] = "    -Box centers for box set subvolume: " + getName();
				lines[3] = "    -Center X, Y, Z, Inline, Crossline";
				lines[4] = "---------------------------------------------------------";
				for (int i = 0; i < getNBoxes(); i++)
				{
					currentBox = (StsBoxSubVolume) (getBoxes().getElement(i));
					center = currentBox.boxCenter;
					double[] xy = project.getAbsoluteXYCoordinates(new float[]
						{center.getX(), center.getY()});
					lines[i + 5] = floatFormat.format(xy[0]) + "   " + floatFormat.format(xy[1]) + "   " + floatFormat.format(center.getZorT()) +
						"   " + currentBox.getRowCenter() + "   " + currentBox.getColCenter();
					if (exportDialog.getAmplitude())
					{
						StsSeismicVolume volume = currentBox.getVolume();
						lines[i + 5] = lines[i + 5] + "   " + volume.getName() + "   " + volume.getIntValue(new float[]
							{center.getX(), center.getY(), center.getZorT()});
					}
				}
				out.writeStringsToFile(lines);
				new StsMessage(currentModel.win3d, StsMessage.INFO, "Successfully output " + filename);
			}
			else
			{
				String[] lines = new String[getNBoxes() + 5];
				lines[0] = "----------------------------------------------------------------------------";
				lines[1] = " Exported from S2S Systems ";
				lines[2] = "    -Boxes for box set subvolume: " + getName();
				lines[3] = "    -Center X, Y, Z, Inline, Crossline, Size X, Y, Z";
				lines[4] = "-----------------------------------------------------------------------------";
				for (int i = 0; i < getNBoxes(); i++)
				{
					currentBox = (StsBoxSubVolume) (getBoxes().getElement(i));
					center = currentBox.boxCenter;
					double[] xy = project.getAbsoluteXYCoordinates(new float[]
						{center.getX(), center.getY()});

					lines[i + 5] = floatFormat.format(xy[0]) + "   " + floatFormat.format(xy[1]) + "   " + floatFormat.format(center.getZorT()) +
						"   " + currentBox.getRowCenter() + "   " + currentBox.getColCenter() + "   " +
						floatFormat.format( (float) ( (float) currentBox.nRows * currentBox.yInc)) + "   " +
						floatFormat.format( (float) ( (float) currentBox.nCols * currentBox.xInc)) + "   " +
						floatFormat.format( (float) ( (float) currentBox.nSlices * currentBox.zInc));
				}
				out.writeStringsToFile(lines);
				new StsMessage(currentModel.win3d, StsMessage.WARNING, "Successfully output " + filename);
			}
		}
		// XML Output
		else
		{
			if (exportDialog.getScope() == exportDialog.CENTERS)
			{
				StsPoint point = null;
				for (int i = 0; i < getNBoxes(); i++)
				{
					currentBox = (StsBoxSubVolume) (getBoxes().getElement(i));
					center = currentBox.boxCenter;
					double[] xy = project.getAbsoluteXYCoordinates(new float[]
						{center.getX(), center.getY()});
					point = new StsPoint(xy[0], xy[1], center.getZorT());
					if (i == 0)
					{
						StsToolkit.writeObjectXML(point, filename, false);
					}
					else
					{
						StsToolkit.writeObjectXML(point, filename, true);
					}
				}
				new StsMessage(currentModel.win3d, StsMessage.INFO, "Successfully output box centers file, " + filename);
			}
			else
			{
				new StsMessage(currentModel.win3d, StsMessage.WARNING, "Output of XML formatted boxes is not currently supportted");
				return false;

				/*                StsOutputBox outBox = new StsOutputBox();
								for(int i=0; i<getNBoxes(); i++)
								{
									currentBox = (StsBoxSubVolume)(getBoxes().getElement(i));
									center = currentBox.boxCenter;
									double[] xy = project.getAbsoluteXYCoordinates(new float[] {center.getX(), center.getY()});
									outBox.center = new StsPoint(xy[0],xy[1],center.getZ());
									outBox.xDim = currentBox.nCols * currentBox.xInc;
									outBox.yDim = currentBox.nRows * currentBox.yInc;
									outBox.zDim = currentBox.nCroppedSlices * currentBox.zInc;
									if(i == 0)
										StsToolkit.writeObjectXML(outBox, filename, false);
									else
										StsToolkit.writeObjectXML(outBox, filename, true);
								}
								new StsMessage(currentModel.win3d, StsMessage.INFO, "Successfully output boxes file, " + filename);
				 */
			}

		}
		return true;
	}

}
