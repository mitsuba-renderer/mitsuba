#include <mitsuba/mitsuba.h>
#include <ply/ply_parser.hpp>

bool ply::ply_parser::parse(std::istream& istream)
{
  std::string line;
  line_number_ = 0;

  std::size_t number_of_format_statements = 0, number_of_element_statements = 0, number_of_property_statements = 0, number_of_obj_info_statements = 0, number_of_comment_statements = 0;

  format_type format = ascii_format;
  std::vector< std::shared_ptr<element> > elements;

  // magic
  char magic[3];
  istream.read(magic, 3);
  istream.ignore(1);
  ++line_number_;
  if (!istream) {
    if (error_callback_) {
      error_callback_(line_number_, "parse error");
    }
    return false;
  }
  if ((magic[0] != 'p') || (magic[1] != 'l') || (magic[2] != 'y')) {
    if (error_callback_) {
      error_callback_(line_number_, "parse error");
    }
    return false;
  }
  if (magic_callback_) {
    magic_callback_();
  }

  while (std::getline(istream, line)) {
    ++line_number_;
    std::istringstream stringstream(line);
    stringstream.unsetf(std::ios_base::skipws);

    stringstream >> std::ws;
    if (stringstream.eof()) {
      if (warning_callback_) {
        warning_callback_(line_number_, "ignoring line '" + line + "'");
      }
    }
    else {
      std::string keyword;
      stringstream >> keyword;

      // format
      if (keyword == "format") {
        std::string format_string, version;
        char space_format_format_string, space_format_string_version;
        stringstream >> space_format_format_string >> std::ws >> format_string >> space_format_string_version >> std::ws >> version >> std::ws;
        if (!stringstream.eof() || !std::isspace(space_format_format_string) || !std::isspace(space_format_string_version)) {
          if (error_callback_) {
            error_callback_(line_number_, "parse error");
          }
          return false;
        }
        if (format_string == "ascii") {
          format = ascii_format;
        }
        else if (format_string == "binary_big_endian") {
          format = binary_big_endian_format;
        }
        else if (format_string == "binary_little_endian") {
          format = binary_little_endian_format;
        }
        else {
          if (error_callback_) {
            error_callback_(line_number_, "parse error");
          }
          return false;
        }
        if (version != "1.0") {
          if (error_callback_) {
            error_callback_(line_number_, "version '" + version + "' is not supported");
          }
          return false;
        }
        if (number_of_format_statements > 0) {
          if (error_callback_) {
            error_callback_(line_number_, "parse error");
          }
          return false;
        }
        ++number_of_format_statements;
        if (format_callback_) {
          format_callback_(format, version);
        }
      }

      // element
      else if (keyword == "element") {
        std::string name;
        std::size_t count;
        char space_element_name, space_name_count;
        stringstream >> space_element_name >> std::ws >> name >> space_name_count >> std::ws >> count >> std::ws;
        if ( !stringstream.eof() || !std::isspace(space_element_name) || !std::isspace(space_name_count)) {
          if (error_callback_) {
            error_callback_(line_number_, "parse error");
          }
          return false;
        }
        std::vector< std::shared_ptr<element> >::const_iterator iterator;
        for (iterator = elements.begin(); iterator != elements.end(); ++iterator) {
          const struct element& element = *(iterator->get());
          if (element.name == name) {
            break;
          }
        }
        if (iterator != elements.end()) {
          if (error_callback_) {
            error_callback_(line_number_, "parse error");
          }
          return false;
        }
        ++number_of_element_statements;
        element_callbacks_type element_callbacks;
        if (element_definition_callbacks_) {
          element_callbacks = element_definition_callbacks_(name, count);
        }
        std::shared_ptr<element> element_ptr(new element(name, count, std::get<0>(element_callbacks), std::get<1>(element_callbacks)));
        elements.push_back(std::shared_ptr<element>(element_ptr));
        current_element_ = element_ptr.get();
      }

      // property
      else if (keyword == "property") {
        std::string type_or_list;
        char space_property_type_or_list;
        stringstream >> space_property_type_or_list >> std::ws >> type_or_list;
        if ( !std::isspace(space_property_type_or_list)) {
          if (error_callback_) {
            error_callback_(line_number_, "parse error");
          }
          return false;
        }
        if (type_or_list != "list") {
          std::string name;
          std::string& type = type_or_list;
          char space_type_name;
          stringstream >> space_type_name >> std::ws >> name >> std::ws;
          if (!stringstream.eof() || !std::isspace(space_type_name)) {
            if (error_callback_) {
              error_callback_(line_number_, "parse error");
            }
            return false;
          }
          if (number_of_element_statements == 0) {
            if (error_callback_) {
              error_callback_(line_number_, "parse error");
            }
            return false;
          }
          std::vector< std::shared_ptr<property> >::const_iterator iterator;
          for (iterator = current_element_->properties.begin(); iterator != current_element_->properties.end(); ++iterator) {
            const struct property& property = *(iterator->get());
            if (property.name == name) {
              break;
            }
          }
          if (iterator != current_element_->properties.end()) {
            if (error_callback_) {
              error_callback_(line_number_, "parse error");
            }
            return false;
          }
          if ((type == type_traits<int8>::name()) || (type == type_traits<int8>::old_name())) {
            parse_scalar_property_definition<int8>(name);
          }
          else if ((type == type_traits<int16>::name()) || (type == type_traits<int16>::old_name())) {
            parse_scalar_property_definition<int16>(name);
          }
          else if ((type == type_traits<int32>::name()) || (type == type_traits<int32>::old_name())) {
            parse_scalar_property_definition<int32>(name);
          }
          else if ((type == type_traits<uint8>::name()) || (type == type_traits<uint8>::old_name())) {
            parse_scalar_property_definition<uint8>(name);
          }
          else if ((type == type_traits<uint16>::name()) || (type == type_traits<uint16>::old_name())) {
            parse_scalar_property_definition<uint16>(name);
          }
          else if ((type == type_traits<uint32>::name()) || (type == type_traits<uint32>::old_name())) {
            parse_scalar_property_definition<uint32>(name);
          }
          else if ((type == type_traits<float32>::name()) || (type == type_traits<float32>::old_name())) {
            parse_scalar_property_definition<float32>(name);
          }
          else if ((type == type_traits<float64>::name()) || (type == type_traits<float64>::old_name())) {
            parse_scalar_property_definition<float64>(name);
          }
          else {
            if (error_callback_) {
              error_callback_(line_number_, "parse error");
            }
            return false;
          }
          ++number_of_property_statements;
        }
        else {
          std::string name;
          std::string size_type_string, scalar_type_string;
          char space_list_size_type, space_size_type_scalar_type, space_scalar_type_name;
          stringstream >> space_list_size_type >> std::ws >> size_type_string >> space_size_type_scalar_type >> std::ws >> scalar_type_string >> space_scalar_type_name >> std::ws >> name >> std::ws;
          if ( !std::isspace(space_list_size_type) || !std::isspace(space_size_type_scalar_type) || !std::isspace(space_scalar_type_name)) {
            if (error_callback_) {
              error_callback_(line_number_, "parse error");
            }
            return false;
          }
          if (number_of_element_statements == 0) {
            if (error_callback_) {
              error_callback_(line_number_, "parse error");
            }
            return false;
          }
          std::vector< std::shared_ptr<property> >::const_iterator iterator;
          for (iterator = current_element_->properties.begin(); iterator != current_element_->properties.end(); ++iterator) {
            const struct property& property = *(iterator->get());
            if (property.name == name) {
              break;
            }
          }
          if (iterator != current_element_->properties.end()) {
            if (error_callback_) {
              error_callback_(line_number_, "parse error");
            }
            return false;
          }
          if (
                (size_type_string == type_traits<uint8>::name()) || (size_type_string == type_traits<uint8>::old_name()) ||
                (size_type_string == type_traits<int8>::name()) || (size_type_string == type_traits<int8>::old_name())
            ) {
            typedef uint8 size_type;
            if ((scalar_type_string == type_traits<int8>::name()) || (scalar_type_string == type_traits<int8>::old_name())) {
              parse_list_property_definition<size_type, int8>(name);
            }
            else if ((scalar_type_string == type_traits<int16>::name()) || (scalar_type_string == type_traits<int16>::old_name())) {
              parse_list_property_definition<size_type, int16>(name);
            }
            else if ((scalar_type_string == type_traits<int32>::name()) || (scalar_type_string == type_traits<int32>::old_name())) {
              parse_list_property_definition<size_type, int32>(name);
            }
            else if ((scalar_type_string == type_traits<uint8>::name()) || (scalar_type_string == type_traits<uint8>::old_name())) {
              parse_list_property_definition<size_type, uint8>(name);
            }
            else if ((scalar_type_string == type_traits<uint16>::name()) || (scalar_type_string == type_traits<uint16>::old_name())) {
              parse_list_property_definition<size_type, uint16>(name);
            }
            else if ((scalar_type_string == type_traits<uint32>::name()) || (scalar_type_string == type_traits<uint32>::old_name())) {
              parse_list_property_definition<size_type, uint32>(name);
            }
            else if ((scalar_type_string == type_traits<float32>::name()) || (scalar_type_string == type_traits<float32>::old_name())) {
              parse_list_property_definition<size_type, float32>(name);
            }
            else if ((scalar_type_string == type_traits<float64>::name()) || (scalar_type_string == type_traits<float64>::old_name())) {
              parse_list_property_definition<size_type, float64>(name);
            }
            else {
              if (error_callback_) {
                error_callback_(line_number_, "parse error");
              }
              return false;
            }
          }
          else if (
                (size_type_string == type_traits<uint16>::name()) || (size_type_string == type_traits<uint16>::old_name()) ||
                (size_type_string == type_traits<int16>::name()) || (size_type_string == type_traits<int16>::old_name())
            ) {
            typedef uint16 size_type;
            if ((scalar_type_string == type_traits<int8>::name()) || (scalar_type_string == type_traits<int8>::old_name())) {
              parse_list_property_definition<size_type, int8>(name);
            }
            else if ((scalar_type_string == type_traits<int16>::name()) || (scalar_type_string == type_traits<int16>::old_name())) {
              parse_list_property_definition<size_type, int16>(name);
            }
            else if ((scalar_type_string == type_traits<int32>::name()) || (scalar_type_string == type_traits<int32>::old_name())) {
              parse_list_property_definition<size_type, int32>(name);
            }
            else if ((scalar_type_string == type_traits<uint8>::name()) || (scalar_type_string == type_traits<uint8>::old_name())) {
              parse_list_property_definition<size_type, uint8>(name);
            }
            else if ((scalar_type_string == type_traits<uint16>::name()) || (scalar_type_string == type_traits<uint16>::old_name())) {
              parse_list_property_definition<size_type, uint16>(name);
            }
            else if ((scalar_type_string == type_traits<uint32>::name()) || (scalar_type_string == type_traits<uint32>::old_name())) {
              parse_list_property_definition<size_type, uint32>(name);
            }
            else if ((scalar_type_string == type_traits<float32>::name()) || (scalar_type_string == type_traits<float32>::old_name())) {
              parse_list_property_definition<size_type, float32>(name);
            }
            else if ((scalar_type_string == type_traits<float64>::name()) || (scalar_type_string == type_traits<float64>::old_name())) {
              parse_list_property_definition<size_type, float64>(name);
            }
            else {
              if (error_callback_) {
                error_callback_(line_number_, "parse error");
              }
              return false;
            }
          }
          else if (
                (size_type_string == type_traits<uint32>::name()) || (size_type_string == type_traits<uint32>::old_name()) ||
                (size_type_string == type_traits<int32>::name()) || (size_type_string == type_traits<int32>::old_name())
            ) {
            typedef uint32 size_type;
            if ((scalar_type_string == type_traits<int8>::name()) || (scalar_type_string == type_traits<int8>::old_name())) {
              parse_list_property_definition<size_type, int8>(name);
            }
            else if ((scalar_type_string == type_traits<int16>::name()) || (scalar_type_string == type_traits<int16>::old_name())) {
              parse_list_property_definition<size_type, int16>(name);
            }
            else if ((scalar_type_string == type_traits<int32>::name()) || (scalar_type_string == type_traits<int32>::old_name())) {
              parse_list_property_definition<size_type, int32>(name);
            }
            else if ((scalar_type_string == type_traits<uint8>::name()) || (scalar_type_string == type_traits<uint8>::old_name())) {
              parse_list_property_definition<size_type, uint8>(name);
            }
            else if ((scalar_type_string == type_traits<uint16>::name()) || (scalar_type_string == type_traits<uint16>::old_name())) {
              parse_list_property_definition<size_type, uint16>(name);
            }
            else if ((scalar_type_string == type_traits<uint32>::name()) || (scalar_type_string == type_traits<uint32>::old_name())) {
              parse_list_property_definition<size_type, uint32>(name);
            }
            else if ((scalar_type_string == type_traits<float32>::name()) || (scalar_type_string == type_traits<float32>::old_name())) {
              parse_list_property_definition<size_type, float32>(name);
            }
            else if ((scalar_type_string == type_traits<float64>::name()) || (scalar_type_string == type_traits<float64>::old_name())) {
              parse_list_property_definition<size_type, float64>(name);
            }
            else {
              if (error_callback_) {
                error_callback_(line_number_, "parse error");
              }
              return false;
            }
          }
          else {
            if (error_callback_) {
              error_callback_(line_number_, "parse error");
            }
            return false;
          }
          ++number_of_property_statements;
        }
      }

      // comment
      else if (keyword == "comment") {
        if (comment_callback_) {
          comment_callback_(line);
        }
        ++number_of_comment_statements;
      }

      // obj_info
      else if (keyword == "obj_info") {
        if (obj_info_callback_) {
          obj_info_callback_(line);
        }
        ++number_of_obj_info_statements;
      }

      // end_header
      else if (keyword == "end_header") {
        if (end_header_callback_) {
          if (end_header_callback_() == false) {
            return true;
          }
        }
        break;
      }

      // unknown keyword
      else {
        if (warning_callback_) {
          warning_callback_(line_number_, "ignoring line '" + line + "'");
        }
      }

    }
  }

  if (number_of_format_statements == 0) {
    if (error_callback_) {
      error_callback_(line_number_, "parse error");
    }
    return false;
  }

  // ascii
  if (format == ascii_format) {
    for (std::vector< std::shared_ptr<element> >::const_iterator element_iterator = elements.begin(); element_iterator != elements.end(); ++element_iterator) {
      struct element& element = *(element_iterator->get());
      for (std::size_t element_index = 0; element_index < element.count; ++element_index) {
        if (element.begin_element_callback) {
          element.begin_element_callback();
        }
        if (!std::getline(istream, line)) {
          if (error_callback_) {
            error_callback_(line_number_, "parse error");
          }
          return false;
        }
        ++line_number_;
        std::istringstream stringstream(line);
        stringstream.unsetf(std::ios_base::skipws);
        stringstream >> std::ws;
        for (std::vector< std::shared_ptr<property> >::const_iterator property_iterator = element.properties.begin(); property_iterator != element.properties.end(); ++property_iterator) {
          struct property& property = *(property_iterator->get());
          if (property.parse(*this, format, stringstream) == false) {
            return false;
          }
        }
        if (!stringstream.eof()) {
          if (error_callback_) {
            error_callback_(line_number_, "parse error");
          }
          return false;
        }
        if (element.end_element_callback) {
          element.end_element_callback();
        }
      }
    }
    istream >> std::ws;
    if (istream.fail() || !istream.eof() || istream.bad()) {
      if (error_callback_) {
        error_callback_(line_number_, "parse error");
      }
      return false;
    }
    return true;
  }

  // binary
  else {
    for (std::vector< std::shared_ptr<element> >::const_iterator element_iterator = elements.begin(); element_iterator != elements.end(); ++element_iterator) {
      struct element& element = *(element_iterator->get());
      for (std::size_t element_index = 0; element_index < element.count; ++element_index) {
        if (element.begin_element_callback) {
          element.begin_element_callback();
        }
        for (std::vector< std::shared_ptr<property> >::const_iterator property_iterator = element.properties.begin(); property_iterator != element.properties.end(); ++property_iterator) {
          struct property& property = *(property_iterator->get());
          if (property.parse(*this, format, istream) == false) {
            return false;
          }
        }
        if (element.end_element_callback) {
          element.end_element_callback();
        }
      }
    }
    if (istream.fail() || (istream.rdbuf()->sgetc() != std::char_traits<char>::eof()) || istream.bad()) {
      if (error_callback_) {
        error_callback_(line_number_, "parse error");
      }
      return false;
    }
    return true;
  }
}
