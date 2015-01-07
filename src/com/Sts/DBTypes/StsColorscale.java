//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.nio.*;
import java.util.*;

/** StsColorscale is the persisted object which contains all relevant colorscale information.
 *  The corresponding UI for StsColorscale is the StsColorscalePanel which has the button and slider controls.
 *
 *  A view-only StsColorscalePanel can be added to a StsColorscaleFieldBean (not editable) for inclusion on the object panel.
 *  Or an editable StsColorscalePanel can be added to an StsEditableColorscalePanel which
 *  includes an edit button for inclusion on the object Panel.
 *
 *  An editable StsColorscalePanel can be inserted into an StsColorscaleDialog for an editable dialog UI.
 *
 *  Colorscale has 255 colorValues from indices 0-254 and 1 null color value (255).  For a given data range,
 *  the first and last colors are half-increments and intermediate colors are full increments, so the total span
 *  is 254 increments.  The color index is therefore Math.round(254*(value - minValue)/(maxValue - minValue)).
 */

public class StsColorscale extends StsMainObject
{
	/** colors used */
	StsSpectrum spectrum;
	/** number of colors */
	int nColors = 255;
	/** total range (min and max) for colorscale */
	float[] range = new float[2];
	/** currently selected range (typically with sliders on a colorscalePanel) */
	float[] editRange = new float[2];
	/** The opacity values for the current spectrum */
	int[] opacityValues = new int[255];
	/** values off of above and below edit range are isVisible as transparent or end color if false */
	boolean transparencyMode = true;
    /** transparent color is replaced with an actual color */
    boolean transparentColorFilled = false;
    /** the full spectrum of colors is compressed between edit range limits if true */
	int compressionMode = CLIPPED;
	/** colors are flipped */
	boolean flip = false;
	/** sliders are synched */
	boolean synched = false;
	/** colors are rotated around by this index amount */
	int rotateAmount = 0;
	/** Voxel keys */
	VoxelKey[] voxelKeys = null;
	/** initial colors defined */
	transient Color[] originalColors = null;
	/** edited (clipped, compressed, etc) from original */
	transient Color[] newColors = null;
	/** min color index isVisible */
	transient int minIndex = 0;
	/** max non-null color index isVisible */
	transient int maxIndex = 255;
    /** color components for scale */
    transient float[][] arrayRGBA = null;
    /** color components in a buffer for shader */
    transient FloatBuffer colormapBuffer = null;
    /** null color currently used (if any) */
    transient public Color nullColor;

    transient private Vector itemListeners = null;
	transient private Vector actionListeners = null;

	transient private boolean voxelsChanged = false;
	transient private boolean opacityChanged = false;

	static final boolean debug = false;

	// ItemEvent states optionally used by ItemListeners
	public static final int CLIPPED = 0;
	public static final int COMPRESSED = 1;
	public static final int COLORS_CHANGED = 2;
	public static final int SOME_COLORS_CHANGED = 3;
	public static final int VOXELS_CHANGED = 4;

	public StsColorscale()
	{
	}

	public StsColorscale(boolean persistent)
	{
		super(persistent);
	}

	public StsColorscale(String name, StsSpectrum spectrum, float[] range)
	{
		super(false);
		setName(name);
		setRange(range);
		for (int i = 0; i < 255; i++)
			opacityValues[i] = 255;
		setSpectrum(spectrum);
		//voxelKeys = StsObjectRefList.constructor(2, 2, "voxelKeys", this);

		initialize();
		addToModel();
	}

	public StsColorscale(String name, StsSpectrum spectrum, float min, float max)
	{
		this(name, spectrum, new float[] {min, max});
	}

    public StsColorscale(StsSpectrum spectrum)
    {
        this(spectrum.name, spectrum, new float[] {0.0f, 0.0f});
    }

	public boolean initialize(StsModel model)
	{
		return initialize();
	}

	public boolean initialize()
	{
		if (spectrum != null)
		{
			initColors(spectrum);
			adjustSpectrum();
		}
		return true;
	}

	// Set colorscale to same settings as passed in colorscale
	public void copySettingsFrom(StsColorscale cs)
	{
		setSpectrum(cs.getSpectrum());
		setCompressionMode(cs.getCompressionMode());
		setEditRange(cs.getEditRange());
		setFlip(cs.getFlip());
		setIsSynched(cs.getIsSynched());
		setTransparencyMode(cs.getTransparencyMode());
		setRotateAmount(cs.getRotateAmount());
		setOpacityValues(cs.opacityValues);
		setNewColors(cs.getNewColors());
		adjustSpectrum();
	}

	/**
	 * Add a voxel key point
	 */
	public VoxelKey addVoxelKey(int min, int max)
	{
		VoxelKey voxelKey = new VoxelKey(min, max);
		voxelKeys = (VoxelKey[])StsMath.arrayAddElement(voxelKeys, voxelKey);
		voxelsChanged();
		return voxelKey;
	}

	/**
	 * Delete a voxel key point
	 */
	public VoxelKey deleteVoxelKey(VoxelKey key)
	{
		if (voxelKeys == null)return null;
		voxelKeys = (VoxelKey[])StsMath.arrayDeleteElement(voxelKeys, key);
		voxelsChanged();

		if (voxelKeys.length == 0)
			return null;
		else
			return voxelKeys[0];
	}

	/**
	 * Edit a voxel key point
	 */
	public void editVoxelKey(VoxelKey findKey, int minIdx, int maxIdx)
	{
		if (voxelKeys == null)return;

		for (int n = 0; n < voxelKeys.length; n++)
		{
			VoxelKey key = voxelKeys[n];
			if (key == findKey)
			{
				key.setMax(maxIdx);
				key.setMin(minIdx);
				voxelsChanged();
				return;
			}
		}
	}

	private void voxelsChanged()
	{
		voxelsChanged = true;
//		dbFieldChanged("voxelKeys", voxelKeys);
		ActionEvent actionEvent = new ActionEvent(this, 0, "voxelsChanged", VOXELS_CHANGED);
		fireActionPerformed(actionEvent);
//		ItemEvent event = new ItemEvent(new JCheckBox(), VOXELS_CHANGED, this, ItemEvent.ITEM_STATE_CHANGED);
//		fireItemStateChanged(event);
	}

	public boolean canDrawVoxels()
	{
		return voxelKeys != null && voxelKeys.length > 0;
	}

	// Spectrum Methods
	public void setSpectrum(StsSpectrum spectrum)
	{
		if(this.spectrum == spectrum) return;
		this.spectrum = spectrum;
		initColors(spectrum);

		minIndex = 0;
		maxIndex = nColors - 1;
		adjustSpectrum();
	}

	public StsSpectrum getSpectrum()
	{
		return spectrum; }

	// Range Methods
	public boolean setRange(float[] range)
	{
		if (range == null || range.length < 2) return false;
	    return setRange(range[0], range[1]);
	}

	public boolean setRange(float newMin, float newMax)
	{
        if(range[0] == newMin && range[1] == newMax) return false;
        range[0] = newMin;
		range[1] = newMax;
        if(debug) System.out.println("StsColorscale.setRange:Edit Range:" + newMin + ", " + newMax);
		return setEditRange(newMin, newMax);
    }

	public float[] getRange()
	{
		return range;
    }

	public boolean setEditRange(float[] range)
	{
		if (editRange == null || editRange.length < 2) return false;
        return setEditRange(range[0], range[1]);
	}

	public boolean setEditRange(float newMin, float newMax)
	{
		editRange[0] = newMin;
		editRange[1] = newMax;

		boolean changed = setIndices(newMin, newMax);
		if(changed)
		{
			adjustSpectrum();
		}
		return changed;
	}

	public boolean setEditMin(float newMin)
	{
		editRange[0] = newMin;
		int newMinIndex = getIndexFromValue(newMin);
		if(newMinIndex == minIndex) return false;
		minIndex = newMinIndex;
		return true;
	}

	public float getEditMin() { return editRange[0]; }

	public boolean setEditMax(float newMax)
	{
		editRange[1] = newMax;
		int newMaxIndex = getIndexFromValue(newMax);
		if(newMaxIndex == maxIndex) return false;
		maxIndex = newMaxIndex;
		return true;
	}

	public float getEditMax() { return editRange[1]; }

	public boolean setIndices(float min, float max)
	{
		// Calculate the new Index....
		boolean changed = false;
		int newMinIndex = getIndexFromValue(min);
		if(newMinIndex != minIndex)
		{
			changed = true;
			minIndex = newMinIndex;
		}
		int newMaxIndex = getIndexFromValue(max);
		if(newMaxIndex != maxIndex)
		{
			changed = true;
			maxIndex = newMaxIndex;
		}
		return changed;
//		float scale = (range[1] - range[0]) / nColors;
//		minIndex = (int) (nColors - (range[1] - min) / scale);
//		maxIndex = (int) (nColors - 1 - (range[1] - max) / scale);
	}

	public int getIndexFromValue(float value)
	{
		if(value <= range[0]) return 0;
		if(value >= range[1]) return 254;
		return Math.round(254*(value - range[0])/(range[1] - range[0]));
	}

	public StsColor getStsColorFromValue(float value)
	{
        int index = getIndexFromValue(value);
        return getStsColor(index);
	}

    public synchronized void removeItemListener(ItemListener listener)
	{
		if (itemListeners != null && itemListeners.contains(listener))
		{
			Vector v = (Vector)itemListeners.clone();
			v.removeElement(listener);
			itemListeners = v;
		}
	}

	public synchronized void addItemListener(ItemListener listener)
	{
		Vector v = itemListeners == null ? new Vector(2) : (Vector)itemListeners.clone();
		if (!v.contains(listener))
		{
			v.addElement(listener);
			itemListeners = v;
		}
	}

	protected void fireItemStateChanged(ItemEvent e)
	{
		if (itemListeners != null)
		{
			Vector listeners = itemListeners;
			int count = listeners.size();
			for(int i = 0; i < count; i++)
			{
				ItemListener listener = (ItemListener)listeners.elementAt(i);
				if(debug)System.out.println("StsColorscale.fireItemStateChanged() listener: " + listener.toString());
				listener.itemStateChanged(e);
			}
		}
	}

	public void itemStateChanged(ItemEvent e)
	{
		ItemEvent event = null;
		fireItemStateChanged(event);
	}

	public synchronized void removeActionListener(ActionListener listener)
	{
		if (actionListeners != null && actionListeners.contains(listener))
		{
			Vector v = (Vector)actionListeners.clone();
			v.removeElement(listener);
			actionListeners = v;
		}
	}

	public synchronized void addActionListener(ActionListener listener)
	{
		Vector v = actionListeners == null ? new Vector(2) : (Vector)actionListeners.clone();
		if (!v.contains(listener))
		{
			v.addElement(listener);
			actionListeners = v;
		}
	}

	protected void fireActionPerformed(ActionEvent e)
	{
		if (actionListeners != null)
		{
			Vector listeners = actionListeners;
			int count = listeners.size();
			for(int i = 0; i < count; i++)
			{
				ActionListener listener = (ActionListener)listeners.elementAt(i);
				if(debug)System.out.println("StsColorscale.fireActionPerformed() listener: " + listener.toString());
				listener.actionPerformed(e);
			}
		}
	}

	public float[] getEditRange()
	{
		return editRange;
	}

	// Mode Methods
	public void setCompressionMode(int mode)
	{
		if(this.compressionMode == mode) return;
		this.compressionMode = mode;
		adjustSpectrum();
		colorsChanged();
	}

	public int getCompressionMode()
	{
		return compressionMode;
	}

	public void setFlip(boolean flip)
	{
		if(this.flip == flip) return;
		this.flip = flip;
		adjustSpectrum();
		colorsChanged();
	}

	public boolean getFlip()
	{
		return flip;
	}

	public boolean getIsSynched() { return synched; }
	public void setIsSynched(boolean synched)
	{
		if(this.synched == synched) return;
		this.synched = synched;
	}

	public void setRotateAmount(int amount)
	{
		rotateAmount = amount;
		adjustSpectrum();
	}

	public int getRotateAmount()
	{
		return rotateAmount;
	}

	// Transparency Methods
	public void setTransparencyMode(boolean mode)
	{
		this.transparencyMode = mode;
		adjustSpectrum();
	}

	public boolean getTransparencyMode()
	{
		return transparencyMode; }

	// Color Methods
	public void initColors(StsSpectrum spectrum)
	{
		if (spectrum == null) return;
		originalColors = spectrum.getColors();
		if (originalColors == null)return;
		newColors = spectrum.getColors();
		nColors = originalColors.length;
		if(nColors != 255)
			StsException.systemError("StsColorscale.initColors() initialized with spectrum which doesn't have 255 colors (one reserved for null).");
		setIndices(editRange[0], editRange[1]);
	}

	public Color[] getOriginalColors()
	{
		Color[] colors = null;
		int i, idx;

		if (flip)
			colors = flipColors(spectrum.getColors());
		else
			colors = spectrum.getColors();

		if (rotateAmount != 0)
		{
			for (i = 0; i < nColors; i++)
			{
				idx = rotateAmount + i;
				if (idx >= nColors - 1)
					idx = idx - nColors;
				colors[i] = colors[idx];
			}
		}
		return colors;
	}

	public Color[] flipColors(Color[] colors)
	{
		Color tColor = null;

		int nColors = colors.length;
		for (int i = 0; i < nColors / 2; i++)
		{
			tColor = colors[i];
			colors[i] = colors[nColors - i - 1];
			colors[nColors - i - 1] = tColor;
		}
//		ItemEvent event = new ItemEvent(new JCheckBox(), COLORS_CHANGED, this, ItemEvent.ITEM_STATE_CHANGED);
//		fireItemStateChanged(event);
		return colors;
	}

	public void setNewColors(Color[] colors)
	{
        arrayRGBA = null;
		this.newColors = colors;
//		colorsChanged();
	}

	public void colorsChanged()
	{
        arrayRGBA = null;
		// Construct and throw item change event
		ActionEvent actionEvent = new ActionEvent(this, MouseEvent.MOUSE_RELEASED, "colorScaleChanged", COLORS_CHANGED);
		fireActionPerformed(actionEvent);
//		ItemEvent event = new ItemEvent(new JCheckBox(), COLORS_CHANGED, this, ItemEvent.ITEM_STATE_CHANGED);
//		fireItemStateChanged(event);
	}

    public void colorsChanged(int eventID)
	{
        arrayRGBA = null;
		// Construct and throw item change event
        if(debug) System.out.println("StsColorscale.colorsChanged() called with eventID " + eventID);
        ActionEvent actionEvent = new ActionEvent(this, eventID, "colorScaleChanged", COLORS_CHANGED);
		fireActionPerformed(actionEvent);
//		ItemEvent event = new ItemEvent(new JCheckBox(), COLORS_CHANGED, this, ItemEvent.ITEM_STATE_CHANGED);
//		fireItemStateChanged(event);
	}

    public void setSomeNewColors(int startIndex, int num, Color color)
	{
		int endIndex = startIndex + num - 1;
//		System.out.println("Setting new colors from " + startIndex + " to " + endIndex + " color " + color.toString() );
		if (startIndex < 0 || startIndex + num >= newColors.length)return;
		for (int i = startIndex; i < startIndex + num; i++)
			this.newColors[i] = color;
		colorsChanged();
	}

	public Color[] getNewColors()
	{
		if (this.newColors == null)
		{
			if (originalColors == null)
				return null;
			newColors = (Color[])StsMath.arraycopy(originalColors, getNColors());
		}
		return this.newColors;
	}

	public Color[] getNewColorsInclTransparency()
	{
		Color[] colors = getNewColors();
		if (colors == null)return null;
		colors = (Color[])StsMath.arrayAddElement(colors, getSpectrum().getTransparentColor().getColor());
		return colors;
	}

	public int getNColors()
	{
		return nColors; }
/*
	public String getName()
	{
		if (spectrum == null)return super.getName() + "-null-spectrum";
		else return super.getName() + spectrum.getName();
	}
*/
	public boolean getUpdateState()
	{
		return true; }

	public void adjustSpectrum()
	{
		int i;

		if(spectrum == null) return;

	    if(debug) System.out.println("StsColorscale.adjustSpectrum(). Colors being adjusted between " + minIndex + " and " + maxIndex);
		Color[] originalColors = getOriginalColors();
		if (originalColors == null) return;
		if (opacityValues != null)
		{
            for (int j = 0; j < nColors; j++)
                if ( (j >= minIndex) && (j <= maxIndex))
                {
                    if(originalColors[j].getAlpha() != opacityValues[j])
                        originalColors[j] = new Color(originalColors[j].getRed(), originalColors[j].getGreen(), originalColors[j].getBlue(), opacityValues[j]);
                }
                else
				{
					if (transparencyMode)
						originalColors[j] = new Color(originalColors[j].getRed(), originalColors[j].getGreen(),
							originalColors[j].getBlue(), 0);
					else
						originalColors[j] = new Color(originalColors[j].getRed(), originalColors[j].getGreen(),
							originalColors[j].getBlue(), 255);
				}
		}
		Color[] newColors = getNewColors();
//        if(newColors == null)
//            newColors = getOriginalColors();

		// Compressed Spectrum
		if (compressionMode == COMPRESSED)
		{
			// Compress the spectrum
			float multiple = (float) (nColors / (float) (maxIndex - minIndex + 1));
			for (i = minIndex; i <= maxIndex; i++)
				newColors[i] = originalColors[ (int) ( (float) (i - minIndex) * multiple)];

			// Set Spectrum Ends
			for (i = nColors - 1; i > maxIndex; i--)
			{
				newColors[i] = originalColors[nColors - 1];
				Color color = newColors[i];
				if (transparencyMode)
					newColors[i] = new Color(color.getRed(), color.getGreen(), color.getBlue(), 0);
				else
					newColors[i] = new Color(color.getRed(), color.getGreen(), color.getBlue(), 255);
			}
			for (i = 0; i < minIndex; i++)
			{
				newColors[i] = originalColors[0];
				Color color = newColors[i];
				if (transparencyMode)
					newColors[i] = new Color(color.getRed(), color.getGreen(), color.getBlue(), 0);
				else
					newColors[i] = new Color(color.getRed(), color.getGreen(), color.getBlue(), 255);
			}
		}
		// Clipped Spectrum
		else
		{
			// Reset Colors between sliders
			for (i = minIndex; i <= maxIndex; i++)
				newColors[i] = originalColors[i];

			// Set Spectrum Ends
			for (i = nColors - 1; i > maxIndex; i--)
			{
				newColors[i] = originalColors[maxIndex];
				Color color = newColors[i];
				if (transparencyMode)
					newColors[i] = new Color(color.getRed(), color.getGreen(), color.getBlue(), 0);
				else
					newColors[i] = new Color(color.getRed(), color.getGreen(), color.getBlue(), 255);

			}
			for (i = 0; i < minIndex; i++)
			{
				newColors[i] = originalColors[minIndex];
				Color color = newColors[i];
				if (transparencyMode)
					newColors[i] = new Color(color.getRed(), color.getGreen(), color.getBlue(), 0);
				else
					newColors[i] = new Color(color.getRed(), color.getGreen(), color.getBlue(), 255);
			}
		}
		setNewColors(newColors);
		// Not sure why this next line was commented out but without it the seismic data comes up wrong
		// until the colorscale editor is entered and jiggled.
		resetOpacity();
	}

	public void resetOpacity()
	{
		if ( (opacityValues == null) || (newColors == null)) return;
		for (int i = minIndex; i < maxIndex; i++)
			newColors[i] = new Color(newColors[i].getRed() / 255.0f, newColors[i].getGreen() / 255.0f,
									 newColors[i].getBlue() / 255.0f, opacityValues[i]/255.0f);
		setNewColors(newColors);
	}

	public void setOpacityValues(int[] opacity)
	{
		opacityValues = opacity;
		opacityChanged = true;
	}

	public int[] getOpacityValues()
	{
		return opacityValues; }

	public Color getColor(int index)
	{
		if (index >= nColors)
			return StsSpectrum.nullColor;
		if (newColors == null)
			return originalColors[index];
		else
			return newColors[index];
	}

	public StsColor getStsColor(int index)
	{
		return new StsColor(getColor(index));
	}

    public Color getNullColor()
    {
        return nullColor;
    }

    /**
	 * Get the data value associated with a particular index
	 * @param idx - the color index being queried.
	 * @return the value associated with the index.
	 */
	public float getDataValue(int idx)
	{
		float scale = (range[1] - range[0]) / nColors;
		float value = idx * scale + range[0];
		return value;
	}

	public void compressColorscale(int percent)
	{
		setCompressionMode(COMPRESSED);
		float curRange = Math.abs(range[1] - range[0]);
		float change = curRange * (float) (percent / 100.0f);
		if ( (editRange[0] + change) < (editRange[1] - change))
        {
            setEditRange(editRange[0] + change, editRange[1] - change);
            if(debug) System.out.println("StsColorscale.compressColorscale:Edit Range:" + (editRange[0] + change) + ", " + (editRange[1] - change));
        }
	}

	public void uncompressColorscale(int percent)
	{
		setCompressionMode(COMPRESSED);
		float curRange = Math.abs(range[1] - range[0]);
		float change = curRange * (float) (percent / 100.0f);
		float min = editRange[0] - change;
		if (min < range[0])min = range[0];
		float max = editRange[1] + change;
		if (max > range[1])max = range[1];
        if(debug) System.out.println("StsColorscale.uncompressColorscale:Edit Range:" + min + ", " + max);
		setEditRange(min, max);
	}

	public void commitChanges(StsColorscale originalColorscale)
	{
		if(voxelsChanged) dbFieldChanged("voxelKeys", voxelKeys);
		if(opacityChanged) dbFieldChanged("opacityValues", opacityValues);
		if(transparencyMode != originalColorscale.transparencyMode)
			dbFieldChanged("transparencyMode", transparencyMode);
		if(compressionMode != originalColorscale.compressionMode)
			dbFieldChanged("compressionMode", compressionMode);
		if(flip != originalColorscale.flip)
			dbFieldChanged("flip", flip);
		if(synched != originalColorscale.synched)
			dbFieldChanged("synched", synched);
		if(rotateAmount != originalColorscale.rotateAmount)
			dbFieldChanged("rotateAmount", rotateAmount);
		if(editRange[0] != originalColorscale.editRange[0] || editRange[1] != originalColorscale.editRange[1])
			dbFieldChanged("editRange", editRange);
		if(spectrum != originalColorscale.spectrum)
			dbFieldChanged("spectrum", spectrum);
	}

	public float[][] computeRGBAArray(boolean nullsFilled)
	{
		if(arrayRGBA != null) return arrayRGBA;
        Color[] colors = getNewColorsInclTransparency();
        //Color[] colors = getNewColors();
        if (nullsFilled)
			colors[colors.length - 1] = new Color(colors[colors.length - 1].getRed(), colors[colors.length - 1].getGreen(),
												  colors[colors.length - 1].getBlue(), 255);

		if (colors == null)return null;
		int nColors = colors.length;
		float norm = 1.0f / 255.0f;

		arrayRGBA = new float[4][nColors];
		float[] rgba = new float[4];
		for (int n = 0; n < nColors; n++)
		{
			colors[n].getComponents(rgba);
			arrayRGBA[0][n] = rgba[0];
			arrayRGBA[1][n] = rgba[1];
			arrayRGBA[2][n] = rgba[2];
			if (rgba[3] >= 0.0f)
				arrayRGBA[3][n] = rgba[3];
			else
				arrayRGBA[3][n] = 1.0f;
		}
		return arrayRGBA;
	}

    public float[][] getComputeRGBAArray(boolean nullsFilled, Color nullColor)
	{
        if(arrayRGBA != null && !nullColorChanged(nullsFilled, nullColor))
            return arrayRGBA;
        return computeRGBAArray(nullsFilled, nullColor);
    }

    public float[][] computeRGBAArray(boolean nullsFilled, Color nullColor)
    {
        float[] rgba = new float[4];
        if(arrayRGBA == null)
        {
            Color[] colors = getNewColorsInclTransparency();
            int nColors = colors.length;

            arrayRGBA = new float[4][nColors];
            for (int n = 0; n < nColors; n++)
            {
                colors[n].getComponents(rgba);
                arrayRGBA[0][n] = rgba[0];
                arrayRGBA[1][n] = rgba[1];
                arrayRGBA[2][n] = rgba[2];
                if (rgba[3] >= 0.0f)
                    arrayRGBA[3][n] = rgba[3];
                else
                    arrayRGBA[3][n] = 1.0f;
            }
        }
        int nColors = arrayRGBA[0].length;
        if (nullsFilled || nullColor != null)
        {           
            if(nullColor == null)
            {
                for(int n = 0; n < 4; n++)
                    arrayRGBA[n][nColors-1] = arrayRGBA[n][nColors-2];
            }
            else
            {
                nullColor.getComponents(rgba);
                for(int n = 0; n < 4; n++)
                    arrayRGBA[n][nColors-1] = rgba[n];
            }
            transparentColorFilled = true;
        }
        else
        {
			int nC = arrayRGBA[0].length; // jbw nColors defined differently globally
            for(int n = 0; n < 4; n++)
                arrayRGBA[n][nC-1] = 0.0f;
            transparentColorFilled = false;
        }
        this.nullColor = nullColor;
        return arrayRGBA;
	}


    public FloatBuffer getComputeColormapBuffer(boolean nullsFilled, Color nullColor)
    {
        if(colormapBuffer != null && !nullColorChanged(nullsFilled, nullColor))
            return colormapBuffer;
        if(arrayRGBA == null)
            arrayRGBA = getComputeRGBAArray(nullsFilled, nullColor);
        if(colormapBuffer == null)
        colormapBuffer = StsJOGLShader.computeColormapBuffer(arrayRGBA, 256);
        return colormapBuffer;
    }

    public int getNTotalColors() { return nColors + 1; }

    public boolean nullColorChanged(boolean nullsFilled, Color nullColor)
    {
        return transparentColorFilled != nullsFilled || this.nullColor != nullColor;
    }

    public VoxelKey[] getVoxelKeys()
	{
		return voxelKeys; }

	public class VoxelKey implements Serializable
	{
		public int min=0, max=0;
		public VoxelKey()
		{
		}

		public boolean initialize(StsModel model)
		{
			return true;
		}

		// CONSTRUCTOR:
		public VoxelKey(int min, int max)
		{
			this.min = min;
			this.max = max;
		}

		// ACCESSORS
		public void setMin(int min)
		{
			this.min = min;
		}

		public int getMin()
		{
			return min;
		}

		public void setMax(int max)
		{
			this.max = max;
		}

		public int getMax()
		{
			return max;
		}

		public float[] getDataRange()
		{
			float floatMin = getDataValue(min);
			float floatMax = getDataValue(max);
			return new float[] {floatMin, floatMax};
		}

		public float getMinValue()
		{
			return getDataValue(min);
		}

		public float getMaxValue()
		{
			return getDataValue(max);
		}

		public void setMinValue(float value)
		{
			min = getIndexFromValue(value);
		}

		public void setMaxValue(float value)
		{
			max = getIndexFromValue(value);
		}
	}

	public void debugPrint()
	{
		System.out.println("Colorscale state: spectrum " + spectrum.getName() + " totalRange: " + this.range[0] + " " + range[1] + " sliderRange: " + editRange[0] + " " + editRange[1]);
	}
}
