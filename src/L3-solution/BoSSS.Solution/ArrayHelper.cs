/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

//import java.lang.reflect.Array;

using System;

namespace BoSSS.Solution.Utils.Formula.Util {

    /**
     * Helper class for work with arrays.
     * @author udav
     */
    class ArrayHelper {
        //**
        // * Searches the specified array of ints for the specified value using enumerative technique.
        // * @param val the value to be searched for.
        // * @param array the array to be searched.
        // * @return index of the search value, if it is contained in array; otherwise -1.
        // */
        //public static int indexOf(int val, int[] array) {
        //    for (int i = 0; i < array.Length; ++i) {
        //        if (val == array[i]) {
        //            return i;
        //        }
        //    }
        //    return -1;
        //}

        /**
         * Searches the specified array of Objects for the specified value using enumerative technique.
         * @param val the value to be searched for. May be null.
         * @param array the array to be searched. Array can contain null values.
         * @return index of the search value, if it is contained in array; otherwise return -1.
         */
        public static int indexOf(Object val, Object[] array) {
            for (int i = 0; i < array.Length; ++i) {
                if (val == array[i] || val != null && val.Equals(array[i])) {
                    return i;
                }
            }
            return -1;
        }

        /**
         * Checks array of ints for specified size. If array`s size is less or equal specified size,
         * then returns the specified array; otherwise returns new array of specified size
         * with old arrays's content.
         * @param array the array to be checket.
         * @param newCapacity the size to be checked.
         * @return If array's size is less or Equals newCapacity, returns the array;
         * otherwise returns new array of the size Equals newCapacity in which first
         * <c>array.Length</c> elements are copied from the array.
         */
        public static int[] ensureCapacity(int[] array, int newCapacity) {
            return ensureCapacity<int>(array, newCapacity);
        }

        /**
         * Checks array of doubles for specified size. If array`s size is less or equal specified size,
         * then returns the specified array; otherwise returns new array of specified size
         * with old arrays's content.
         * @param array the array to be checket.
         * @param newCapacity the size to be checked.
         * @return If array's size is less or Equals newCapacity, returns the array;
         * otherwise returns new array of the size Equals newCapacity in which first
         * <c>array.Length</c> elements are copied from the array.
         */
        public static double[] ensureCapacity(double[] array, int newCapacity) {
            return ensureCapacity<double>(array, newCapacity);
        }

        /**
         * Checks array for specified size. If array`s size is less or equal specified size,
         * then returns the specified array; otherwise returns new array of specified size
         * with old arrays's content.
         * @param array the array to be checket.
         * @param newCapacity the size to be checked.
         * @return If array's size is less or Equals newCapacity, returns the array;
         * otherwise returns new array of the size Equals newCapacity in which first
         * <c>array.Length</c> elements are copied from the array.
         */
        public static T[] ensureCapacity<T>(T[] array, int newCapacity) {
            if (array.Length >= newCapacity) {
                return array;
            }
            Array.Resize<T>(ref array, newCapacity);
            return array;
        }
    }
}