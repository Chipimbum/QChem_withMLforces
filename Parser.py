import re

ORIGIN_LIST_FILE_INIT = "This is the list of origin files that are found in '{}':\n"
START_OF_ORIENTATION = "\s*Standard Nuclear Orientation\s*"
START_OF_ENERGY = "\s*Total energy in the final basis set =\s*(-?\d+\.\d+)\s*"
START_OF_ATOM_LIST = "\$molecule"
END_OF_ATOM_LIST = "$end"
NEW_LINE = "\n"
ORIENTATION_LINE = "\s+\d+\s+(\w+.*)"
START_OF_GRADIENT = " Gradient of SCF Energy"
END_OF_GRADIENT = " RMS gradient           ="
END_OF_FILE = ''
COL_SEPARATOR = '     '
MINUS = '-'
FORCE_HA_BOHR_TO_EV_ANG = 51.422086190832 #http://greif.geo.berkeley.edu/~driver/conversions.html
ENERGY_HA_TO_EV =  27.211386245988 #https://physics.nist.gov/cgi-bin/cuu/Value?hrev

def _get_force(gradient):
    """
    Gets a gradient, negate it and change the units from a.u. to eV.
    :param value: string of value to convert
    :return: string of new value
    """
    return str(-1 * float(gradient) * FORCE_HA_BOHR_TO_EV_ANG)

def _convert_to_eV(energy):
    """
    Gets a string of energy in a.u. units and return a new string in eV units
    :param energy: string of the energy or force
    :return:
    """
    return str(float(energy) * ENERGY_HA_TO_EV)

class Parser:
    """
    Parse an .out file/s and make one .xyz file out of them.
    """

    def __init__(self, output_file, origin_list_file):
        self._input_file = None
        self._line = None
        self._output_file = output_file
        self._atom_num = None
        self._origin_list_file = origin_list_file
        self._origin_list_file.write(ORIGIN_LIST_FILE_INIT.format(output_file.name))

    def parse_file(self):
        """
        Parses the current input file and add the xyz format of it to the output file.
        """
        while self._search_line(START_OF_ORIENTATION):
            orientations = self._get_orientation()
            energy = self._get_energy()
            gradients = self._get_forces()
            if not orientations or not energy or not gradients \
                    or len(orientations) != self._atom_num or len(gradients) != self._atom_num: return
            self._write_atom_num()
            self._write_energy(energy)
            self._write_orientation_and_gradient(orientations, gradients)

    def _get_forces(self):
        """
        Finds the forces of the current samples
        :return: a list of strings which are the lines of the orientations
        """
        self._search_line(START_OF_GRADIENT)
        if self._end_of_file(): return
        # moving to the start of the matrix
        self._advance(2)
        # making a (num of molecules) * 3 matrix of gradients
        gradients = []
        while not self._line.startswith(END_OF_GRADIENT):
            # takes the first gradient of up to 6 first atoms, removing the index in the start of the line
            first = [_get_force(gradient) for gradient in self._line.split()[1:]]
            self._advance()
            # takes the second gradient of up to 6 first atoms, removing the index in the start of the line
            second = self._line.split()[1:]
            # adding the second gradient the the first
            if len(first) != len(second):
                return
            for gradient in range(len(first)):
                first[gradient] = first[gradient] + COL_SEPARATOR + _get_force(second[gradient])
            #print(first)
            self._advance()
            # takes the third gradient of up to 6 first atoms, removing the index in the start of the line
            third = self._line.split()[1:]
            # adding the second gradient the the first
            if len(first) != len(third): return
            for gradient in range(len(first)):
                first[gradient] = first[gradient] + COL_SEPARATOR + _get_force(third[gradient])
            self._advance(2)
            # adding first second third gradients of x molecules to the atom's gradient list
            gradients += first
        return gradients

    def _write_atom_num(self):
        """
        Writes the number of atoms in the samples to the output file
        :return:
        """
        self._output_file.write(str(self._atom_num) + NEW_LINE)

    def _get_energy(self):
        """
        Finds the energy of the current samples
        :return: a string of the energy
        """
        energy_line = self._search_line(START_OF_ENERGY)
        if energy_line:
            return _convert_to_eV(energy_line.group(1))

    def _write_energy(self, energy):
        """
        Writes the energy of the samples to the output file
        :param energy the energy of the current samples
        """
        self._output_file.write(energy + NEW_LINE)

    def _get_orientation(self):
        """
        Finds the orientations of the current samples
        :return: a list of strings which are the lines of the orientations
        """
        self._search_line(START_OF_ORIENTATION)
        if self._end_of_file(): return
        # moving to the start of the matrix
        self._advance(3)
        # writing each line of the matrix without the index
        orientations = []
        for atom in range(self._atom_num):
            orientation_line = re.match(ORIENTATION_LINE, self._line)
            if not orientation_line:
                return
            orientations.append(re.match(ORIENTATION_LINE, self._line).group(1))
            self._advance()
        return orientations

    def _write_orientation_and_gradient(self, orientation_lines, gradient_lines):
        """
        Write the orientation of the atoms matrix to the output file
        """
        for line in range(len(orientation_lines)):
            self._output_file.write(orientation_lines[line] + COL_SEPARATOR + gradient_lines[line] + NEW_LINE)

    def _search_line(self, string_to_search):
        """
        Moves to the line that has the given string
        :param string_to_search:  string to search for
        :return the match object that if the line was found and None otherwise
        """
        while self._line:
            match = re.match(string_to_search, self._line)
            if match:
                return match
            self._advance()

    def _set_atom_num(self):
        """
        Finds the number of atoms in the samples
        """
        self._search_line(START_OF_ATOM_LIST)
        # get to the start of the atoms matrix
        self._advance(2)
        if self._end_of_file():
            return
        # start counting the lines in the matrix
        lines_between_start_to_end = 0
        while (not self._line.startswith(END_OF_ATOM_LIST)) and self._line != NEW_LINE:
            self._advance()
            lines_between_start_to_end += 1
            if self._end_of_file():
                return
        self._atom_num = lines_between_start_to_end

    def set_new_file(self, input_file):
        """
        Setting a new file to the parser
        :param input_file: name of the new file
        """
        self._input_file = input_file
        self._line = input_file.readline()
        self._set_atom_num()
        if not self._atom_num:
            return
        # writing the name of the file the information was taken from to the origin_list_file
        self._origin_list_file.write(self._input_file.name + NEW_LINE)

    def _advance(self, line_num=1):
        """
        Advance to the next line in the input file
        :param line_num: number of lines to advance
        """
        for line in range(line_num):
            self._line = self._input_file.readline()

    def _end_of_file(self):
        """
        :return: True if all of the file has been read and false otherwise
        """
        return self._line == END_OF_FILE

