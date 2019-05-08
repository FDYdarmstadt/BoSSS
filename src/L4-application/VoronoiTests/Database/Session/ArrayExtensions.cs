using System;

namespace VoronoiTests.Database.Session
{
    static class ArrayMethods
    {
        public static Tout[] ConvertArray<Tin, Tout>(Tin[] array)
        {
            Tout[] newArray = new Tout[array.Length];
            for(int i = 0; i < array.Length; ++i)
            {
                newArray[i] = (Tout)Convert.ChangeType(array[i], typeof(Tout));
            }
            return newArray;
        }
    }
}
