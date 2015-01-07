package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.MVC.*;
import com.Sts.UI.*;

import java.awt.*;

public class StsTypeLibrary extends StsMainObject
{
    protected StsObjectRefList types = null;
    protected StsType currentType = null;

    static public String genericLibraryName = "Generic Library";
    static public String defaultLibraryName = "Lithology Library";

    public StsTypeLibrary()
    {
    }

    public StsTypeLibrary(String name)
    {
        setName(new String(name));
        types = StsObjectRefList.constructor(2, 2, "types", this);
    }

    public boolean initialize(StsModel model)
    {
        return true;
    }

    public void addType(String name, StsColor color)
    {
        StsType type = new StsType(name, color);
        types.add(type);
    }

    public StsType getType(String name, Color color)
    {
        return getType(name, new StsColor(color));
    }

    public StsType getType(String name, StsColor color)
    {
        StsType type;

        int nTypes = types.getSize();
        // look in types in this typeLibrary for this name; if found, return
        for (int n = 0; n < nTypes; n++)
        {
            type = (StsType) types.getElement(n);
            if (name == null)
                return type;
            if (type.getName().equals(name))
                return type;
        }
        // look for type with same name in all libraries: if found, add to this library and return
        type = (StsType) currentModel.getObjectWithName(StsType.class, name);
        // if not found, create type, add to library, and return
        if (type == null)
            type = new StsType(name, new StsColor(color));
        types.add(type);
        return type;
    }

    public StsObjectRefList getTypes()
    {
        return types;
    }

    public StsType getType(int index)
    {
        return (StsType) types.getElement(index);
    }

    public String getTypeName(int index)
    {
        StsType type = (StsType) types.getElement(index);
        if (type == null)
            return "none";
        return type.getName();
    }

    public int getIndexOfType(StsType type)
    {
        return types.getIndex(type);
    }

    public StsColor[] getStsColors()
    {
        if (types == null)
            return new StsColor[0];
        int nColors = types.getSize();
        StsColor[] colors = new StsColor[nColors];
        for (int n = 0; n < nColors; n++) {
            StsType type = (StsType) types.getElement(n);
            colors[n] = type.getStsColor();
        }
        return colors;
    }

    public StsColorListItem getTypeItem() { return new StsColorListItem(currentType); }
    public void setTypeItem(StsColorListItem colorListItem)
    {
        currentType = (StsType)colorListItem.getObject();
    }

	public void setCurrentType(int index)
	{
		if(index >= types.getSize())
			index = 0;
		currentType = (StsType)types.getElement(index);
	}

    public StsColorListItem[] getColorListItems()
    {
        int nTypes = types.getSize();
        StsColorListItem[] colorListItems = new StsColorListItem[nTypes];
        for(int n = 0; n < nTypes; n++)
        {
            StsType type = (StsType)types.getElement(n);
            colorListItems[n] = new StsColorListItem(type);
        }
        return colorListItems;
    }

    public StsType getCurrentType() { return currentType; }

    static public StsTypeLibrary getCreate(String name)
    {
        StsTypeLibrary library = (StsTypeLibrary) currentModel.getObjectWithName(StsTypeLibrary.class, name);
        if (library == null) library = new StsTypeLibrary(name);
        return library;
    }

    static public StsTypeLibrary getCreateGenericLibrary()
    {
        StsTypeLibrary library = (StsTypeLibrary) currentModel.getObjectWithName(StsTypeLibrary.class, genericLibraryName);
        if (library == null) library = createGenericLibrary();
        return library;
    }

    static public StsTypeLibrary getCreateDefaultLibrary()
    {
        StsTypeLibrary library = (StsTypeLibrary) currentModel.getObjectWithName(StsTypeLibrary.class, defaultLibraryName);
        if (library == null) library = createDefaultLibrary();
        return library;
    }

    static public StsTypeLibrary createDefaultLibrary()
    {
        StsTypeLibrary typeLibrary = new StsTypeLibrary(defaultLibraryName);
        StsSpectrum spectrum = currentModel.getSpectrum("Basic");
        StsColor[] colors = spectrum.getStsColors();

        typeLibrary.addType("Sandstone", colors[0]);
        typeLibrary.addType("Shale", colors[1]);
        typeLibrary.addType("Carbonate", colors[2]);
        typeLibrary.addType("Oil Sand", colors[3]);
        typeLibrary.addType("Gas Sand", colors[4]);
        typeLibrary.addType("Cap Rock", colors[5]);
        typeLibrary.addType("Water Sand", colors[6]);
        typeLibrary.addType("Channel Sand", colors[7]);
        typeLibrary.addType("Unknown", colors[8]);

		typeLibrary.setCurrentType(0);
        currentModel.instanceChange(typeLibrary, "constructedDefaultLibrary");

        return typeLibrary;
    }

    static public StsTypeLibrary createGenericLibrary()
    {
        StsType type;

        StsTypeLibrary typeLibrary = new StsTypeLibrary(genericLibraryName);
        StsSpectrum spectrum = currentModel.getSpectrum("Basic");
        StsColor[] colors = spectrum.getStsColors();

        typeLibrary.addType("Red", colors[0]);
        typeLibrary.addType("Yellow", colors[1]);
        typeLibrary.addType("Green", colors[2]);
        typeLibrary.addType("Cyan", colors[3]);
        typeLibrary.addType("Blue", colors[4]);
        typeLibrary.addType("Magenta", colors[5]);
        typeLibrary.addType("Orange", colors[6]);
        typeLibrary.addType("Lime", colors[7]);
        typeLibrary.addType("Sea Green", colors[8]);
        typeLibrary.addType("Slate Blue", colors[9]);
        typeLibrary.addType("Purple", colors[10]);
        typeLibrary.addType("Fushia", colors[11]);
        typeLibrary.addType("Dark Red", colors[12]);
        typeLibrary.addType("Khaki", colors[13]);
        typeLibrary.addType("S2S Green", colors[14]);

	    typeLibrary.setCurrentType(0);
        currentModel.instanceChange(typeLibrary, "constructedGenericLibrary");

        return typeLibrary;
    }
}