/**
 * Professional Button Component - IBM-Inspired Design
 *
 * Features:
 * - Sharp rectangular design (no rounded corners)
 * - IBM-inspired color palette
 * - Multiple variants (primary, secondary, danger, ghost)
 * - Multiple sizes (sm, md, lg)
 * - Accessible with focus states
 * - Loading state support
 */

import React from 'react';

export interface ButtonProps extends React.ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: 'primary' | 'secondary' | 'danger' | 'ghost';
  size?: 'sm' | 'md' | 'lg';
  loading?: boolean;
  fullWidth?: boolean;
  children: React.ReactNode;
}

const Button = React.forwardRef<HTMLButtonElement, ButtonProps>(
  (
    {
      variant = 'primary',
      size = 'md',
      loading = false,
      fullWidth = false,
      disabled,
      className = '',
      children,
      ...props
    },
    ref
  ) => {
    // Base classes - sharp corners, no border radius
    const baseClasses = 'sharp font-medium transition-fast focus-visible:outline-2 focus-visible:outline-offset-2';

    // Variant classes
    const variantClasses = {
      primary: 'btn-primary',
      secondary: 'btn-secondary',
      danger: 'btn-danger',
      ghost: 'btn-ghost',
    };

    // Size classes
    const sizeClasses = {
      sm: 'btn-sm',
      md: '', // Default height from btn-primary/secondary
      lg: 'btn-lg',
    };

    // Width classes
    const widthClasses = fullWidth ? 'w-full' : '';

    // Combine all classes
    const buttonClasses = `
      ${baseClasses}
      ${variantClasses[variant]}
      ${sizeClasses[size]}
      ${widthClasses}
      ${className}
    `.trim().replace(/\s+/g, ' ');

    return (
      <button
        ref={ref}
        className={buttonClasses}
        disabled={disabled || loading}
        {...props}
      >
        {loading ? (
          <span className="flex items-center justify-center gap-2">
            <span className="spinner" />
            <span>{children}</span>
          </span>
        ) : (
          children
        )}
      </button>
    );
  }
);

Button.displayName = 'Button';

export default Button;
